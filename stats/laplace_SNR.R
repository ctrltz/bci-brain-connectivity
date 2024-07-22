###
# Statistical analysis of SNR measured at C3- and C4-Laplace
#
# Outputs:
# 1. laplace_SNR.RData - subject-average and original SNR values, mixed models
#    (accuracy ~ SNR, SNR ~ session), correlation between SNR and accuracy,
#    t-test for group differences in SNR
# 2. laplace_SNR.tex - all the statistical results from above exported to TeX
###


# Prefix for TeX output
prefix = "laplace"

# Load the data
results <- load_C3C4Laplace(input_filename)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)



### Between-subject effect of SNR on accuracy

# Extract the average SNR over both channels and all sessions
SNR_rest_avg <- results %>% aggregate(. ~ Subject, mean)
SNR_rest_avg <- merge(SNR_rest_avg, accuracy_avg, by = 'Subject')

# Check the outliers
assert("Checking for outliers in average SNR", 
       sum(checkOutliers(SNR_rest_avg$LogSNR)) == 0)

# Plot the histogram
p_hist <- ggplot(data = SNR_rest_avg, mapping = aes(x = LogSNR)) +
  geom_histogram(color = "white", fill = "grey") + 
  theme_classic()

# Check the difference between groups
SNR_group_diff_ttest <- with(SNR_rest_avg, t.test(LogSNR[Group == 'MBSR'],
                                                  LogSNR[Group == 'Control']))
SNR_group_diff_cohens_d <- cohens_d(SNR_rest_avg, LogSNR ~ Group, var.equal = F)
print(SNR_group_diff_ttest)
print(SNR_group_diff_cohens_d)

# Correlate
(c_SNR_accuracy <- with(SNR_rest_avg, cor.test(LogSNR, Accuracy)))

# Plot
p1 <- ggplot(data = SNR_rest_avg, mapping = aes(x = LogSNR, y = Accuracy)) +
  geom_point() +
  geom_smooth(method = 'lm', se = T, color = "blue") +
  annotate("text", x = 11.5, y = 0.225, 
           label = paste0('ρ = ', round(c_SNR_accuracy$estimate, 2), '\n', 
                          'p ', format.pval(as.numeric(c_SNR_accuracy$p.value), 0.001, digits = 1)),
           hjust = 1, vjust = 0, size = 3.5, fontface = mark(c_SNR_accuracy$p.value)) +
  xlab('SNR (dB)') + ylab('Accuracy') +
  ylim(0.2, 1) + 
  theme_classic()

# Raincloud plot - compare SNR between groups
p_group_cmp <- ggplot(data = SNR_rest_avg,
                      aes(Group, LogSNR, fill = Group)) +
  geom_rain(alpha = .5) +
  geom_signif(comparisons = list(c("Control", "MBSR")), y_position = 12,
              tip_length = 0, test = "t.test", textsize = 4,
              map_signif_level = c('*' = 0.05, 'n.s.' = 1)) +
  theme_classic() +
  scale_fill_brewer(palette = 'Dark2') +
  scale_y_continuous(breaks = seq(0, 12, by = 4)) +
  expand_limits(y = 14) +
  ylab('SNR (dB)') +
  guides(fill = 'none', color = 'none')

# Export
tex.save(output_filename, "% Between-subject effect of SNR on accuracy\n")
tex.save(output_filename, "SubjectAverageSNRrestAccuracyCor", 
         format(as.numeric(c_SNR_accuracy$estimate), digits = 2), prefix = prefix)
tex.save(output_filename, "SubjectAverageSNRrestAccuracyDf", 
         format(as.numeric(c_SNR_accuracy$parameter)), prefix = prefix)
tex.save(output_filename, "SubjectAverageSNRrestAccuracyTvalue", 
         format(as.numeric(c_SNR_accuracy$statistic), digits = 2), prefix = prefix)
tex.save(output_filename, "SubjectAverageSNRrestAccuracyPvalue", 
         format.pval(as.numeric(c_SNR_accuracy$p.value), 0.001, digits = 1), prefix = prefix)
tex.save(output_filename, "SubjectAverageSNRrestAccuracyCIMin", 
         format(as.numeric(c_SNR_accuracy$conf.int[1]), digits = 2), prefix = prefix)
tex.save(output_filename, "SubjectAverageSNRrestAccuracyCIMax", 
         format(as.numeric(c_SNR_accuracy$conf.int[2]), digits = 2), prefix = prefix)

tex.save(output_filename, "\n% Group Differences in Mean SNR\n")
tex.save(output_filename, "GroupDiffSNRDf", 
         format(round(SNR_group_diff_ttest$parameter, 1), nsmall = 1), prefix = prefix)
tex.save(output_filename, "GroupDiffSNRTvalue", 
         format(round(SNR_group_diff_ttest$statistic, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "GroupDiffSNRPvalue", 
         format.pval(round(SNR_group_diff_ttest$p.value, 2), 0.001, digits = 2), prefix = prefix)
tex.save(output_filename, "GroupDiffSNRCIMin", 
         format(round(SNR_group_diff_ttest$conf.int[[1]], 2), digits = 2), prefix = prefix)
tex.save(output_filename, "GroupDiffSNRCIMax", 
         format(round(SNR_group_diff_ttest$conf.int[[2]], 2), digits = 2), prefix = prefix)
# NOTE: order of groups in cohen's d and t-test does not match
tex.save(output_filename, "GroupDiffSNRCohensD", 
         format(round(as.matrix(-SNR_group_diff_cohens_d[['effsize']])[[1,1]], 2), 
                digits = 2), prefix = prefix)


### Within-subject effect of SNR on accuracy

# Extract the SNR during resting-state, average over C3-Lap and C4-Lap, and scale
SNR_rest <- results %>%
  aggregate(. ~ Subject + Session, mean)
SNR_rest$LogSNR_orig <- SNR_rest$LogSNR
SNR_rest$Accuracy_orig <- SNR_rest$Accuracy
SNR_rest$LogSNR <- scale(SNR_rest$LogSNR)
SNR_rest$Accuracy <- scale(SNR_rest$Accuracy)
SNR_rest$Session_orig <- SNR_rest$Session
SNR_rest$Session <- scale(SNR_rest$Session)

# Check the outliers
assert("Checking for outliers in SNR", 
       sum(checkOutliers(SNR_rest$LogSNR)) == 0)

# Plot average SNR for different sessions
p_ses_avg <- SNR_rest %>%
  aggregate(LogSNR_orig ~ Session_orig, 
            FUN = function(x) c(mn = mean(x), se = sd(x) / sqrt(length(x)))) %>%
  do.call(data.frame, .) %>%
  ggplot(data = ., mapping = aes(x = Session_orig, y = LogSNR_orig.mn)) +
  geom_point() +
  geom_pointrange(aes(ymin = LogSNR_orig.mn - 1.96 * LogSNR_orig.se,
                      ymax = LogSNR_orig.mn + 1.96 * LogSNR_orig.se)) +
  geom_smooth(method = 'lm', se = F, color = "blue") +
  scale_x_continuous(breaks = seq(1, 11, by = 2)) +
  xlab('Session') + ylab('SNR (dB)') +
  theme_classic()


# Fit the mixed model
fm_within <- lmer(formula = snr_vs_accuracy_formula,
                  data = SNR_rest)
assert("Accuracy ~ Laplace SNR did not converge",
       has_converged(fm_within))
coef_snr <- lmerTest:::get_coefmat(fm_within)["LogSNR",]
coef_ci_within <- confint(fm_within)
summary(fm_within)

# Restore original betas
beta0 <- lmerTest:::get_coefmat(fm_within)["(Intercept)","Estimate"]
beta1 <- lmerTest:::get_coefmat(fm_within)["LogSNR","Estimate"]
c(beta0_orig, beta1_orig) %<-% restoreBetas(beta0, beta1, SNR_rest, 
                                            x = 'LogSNR_orig', y = 'Accuracy_orig')

# Plot
p2 <- ggplot(data = SNR_rest, mapping = aes(x = LogSNR_orig, y = Accuracy_orig)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(intercept = beta0_orig, slope = beta1_orig, 
              color = "blue", linewidth = 1) +
  annotate("text", x = 15, y = 0.225, 
           label = paste0('β = ', round(beta1, 2), '\n', 
                          'p ', format.pval(as.numeric(coef_snr[["Pr(>|t|)"]]), 0.001, digits = 1)),
           hjust = 1, vjust = 0, size = 3.5, fontface = mark(coef_snr[["Pr(>|t|)"]])) +
  scale_color_grey() +
  xlab('SNR (dB)') + ylab('Accuracy') +
  expand_limits(y = 0.2) + 
  theme_classic()

# Export
tex.save(output_filename, "\n% Within-subject effect of SNR on accuracy\n")
tex.save(output_filename, "WithinSubjectSNRrestAccuracyBeta", 
         format(as.numeric(coef_snr[["Estimate"]]), digits = 2), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestAccuracyDf", 
         format(as.numeric(coef_snr[["df"]]), digits = 3), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestAccuracyTvalue", 
         format(as.numeric(coef_snr[["t value"]]), digits = 3), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestAccuracyPvalue", 
         format.pval(as.numeric(coef_snr[["Pr(>|t|)"]]), 0.001, digits = 1), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestAccuracyBetaCIMin", 
         format(round(coef_ci_within["LogSNR", 1], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestAccuracyBetaCIMax", 
         format(round(coef_ci_within["LogSNR", 2], 2), nsmall = 2), prefix = prefix)


### Longitudinal changes in SNR

# Fit the mixed model
fm_session <- lmer(session_vs_snr_formula, 
                   data = SNR_rest)
assert("Laplace SNR ~ Session did not converge",
       has_converged(fm_session))
coef_session <- lmerTest:::get_coefmat(fm_session)["Session",]
coef_ci_session <- confint(fm_session)
summary(fm_session)

# Restore the original betas
beta0 <- lmerTest:::get_coefmat(fm_session)["(Intercept)","Estimate"]
beta1 <- lmerTest:::get_coefmat(fm_session)["Session","Estimate"]
c(beta0_orig, beta1_orig) %<-% restoreBetas(beta0, beta1, SNR_rest, 
                                            x = 'Session_orig', y = 'LogSNR_orig')

# Plot the result
p3 <- ggplot(data = SNR_rest, mapping = aes(x = Session_orig, y = LogSNR_orig)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(intercept = beta0_orig, slope = beta1_orig, 
              color = 'blue', linewidth = 1) +
  annotate("text", x = 11, y = 15, 
           label = paste0('β = ', round(beta1, 2), '\n', 
                          'p ', format.pval(as.numeric(coef_session[["Pr(>|t|)"]]), 0.001, digits = 2)),
           hjust = 1, vjust = 1, size = 3.5, fontface = mark(coef_session[["Pr(>|t|)"]])) +
  scale_color_grey() +
  scale_x_continuous(breaks = seq(1, 11, by = 2)) +
  xlab('Session') + ylab('SNR (dB)') +
  theme_classic()


# Export
tex.save(output_filename, "\n% Longitudinal Changes in SNR\n")
tex.save(output_filename, "WithinSubjectSNRrestSessionBeta", 
         format(as.numeric(coef_session[["Estimate"]]), digits = 1), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestSessionDf", 
         format(as.numeric(coef_session[["df"]]), digits = 3), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestSessionTvalue", 
         format(as.numeric(coef_session[["t value"]]), digits = 3), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestSessionPvalue", 
         format.pval(as.numeric(coef_session[["Pr(>|t|)"]]), 0.001, digits = 2), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestSessionBetaCIMin", 
         format(round(coef_ci_session["Session", 1], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "WithinSubjectSNRrestSessionBetaCIMax", 
         format(round(coef_ci_session["Session", 2], 2), nsmall = 2), prefix = prefix)


### Plot examples of different SNR
example_data <- sanity$example.snr[,,1]
subj_low_snr <- example_data$subj.low.snr
subj_med_snr <- example_data$subj.med.snr
subj_high_snr <- example_data$subj.high.snr
spec.avg <- 10 * log10(example_data$spec.avg)
snr.avg <- example_data$snr.avg  # log-scale is applied in MATLAB before averaging
example_df <- data.frame(freqs = example_data$freqs, lowSNR = spec.avg[subj_low_snr,],
                         medSNR = spec.avg[subj_med_snr,], highSNR = spec.avg[subj_high_snr,])
colors <- c("LowSNR" = '#fcbba1', "MedSNR" = '#ef3b2c', "HighSNR" = '#67000d')
labels <- c("LowSNR" = sprintf("SNR = %.2f dB", snr.avg[subj_low_snr,]),
            "MedSNR" = sprintf("SNR = %.2f dB", snr.avg[subj_med_snr,]),
            "HighSNR" = sprintf("SNR = %.2f dB", snr.avg[subj_high_snr,]))
p_examples <- ggplot(example_df, aes(x = freqs)) +
  geom_line(aes(y = lowSNR, color = "LowSNR"), linewidth = 1) + 
  geom_line(aes(y = medSNR, color = "MedSNR"), linewidth = 1) + 
  geom_line(aes(y = highSNR, color = "HighSNR"), linewidth = 1) + 
  scale_colour_manual(values = colors, labels = labels,
                      limits=c("HighSNR", "MedSNR", "LowSNR")) +
  xlim(1, 40) + ylim(-25, 5) +
  xlab('Frequency (Hz)') + ylab(expression('10'%.%'log'[10]*'(PSD)')) +
  theme_classic() + 
  theme(legend.position = c(0.7, 0.8),
        legend.title = element_blank())


### Save the plots together to align them later with the multiverse ones
p_laplace <- list(between = p1, within = p2, session = p3,
                  group_cmp = p_group_cmp, dynamics = p_ses_avg, examples = p_examples)



### Plot individual relationships between SNR and accuracy
p_individual <- ggplot(data = SNR_rest, mapping = aes(x = LogSNR, y = Accuracy_orig)) +
  geom_point(size = 1) + geom_smooth(method = "lm") +
  facet_wrap(. ~ Subject, nrow = 8) +
  expand_limits(y = 0.2) +
  xlab("SNR (dB)") + ylab("Accuracy") + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
ggsave(file.path(plot.path, 
                 paste(prefix, '_snr_rest_vs_accuracy_subject_wise.png', sep = '')),
       plot = p_individual, width = 10, height = 8)


### Save all the results
save(SNR_rest_avg, SNR_rest, SNR_group_diff_ttest, SNR_group_diff_cohens_d,
     c_SNR_accuracy, fm_within, coef_ci_within, fm_session, coef_ci_session,
     file = file.path(r.path, "laplace_SNR.RData"))