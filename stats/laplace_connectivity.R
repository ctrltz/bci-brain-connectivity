###
# Statistical analysis of PS between C3-, C4-, CP3-, and CP4-Laplace channels
#
# Outputs:
# 1. laplace_connectivity.RData - within-subject effects of SNR on PS, PS on
#    BCI performance, and longitudinal changes in PS
# 2. figCsupp1-laplace-snr-connectivity.png - effect of SNR on PS 
# 3. figCsupp2-laplace-connectivity-performance-within - within-subject effect
#      of PS on BCI performance
# 4. figCsupp3-laplace-connectivity-longitude.png - longitudinal changes in PS
# 5. laplace_connectivity.tex - all the statistical results from above 
#    exported to TeX
###

# Create output folder for images
prefix = 'laplace_connectivity'
output.folder <- paste('../results/stats/', prefix, sep='')
dir.create(output.folder)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)

# Load the data
data <- readMat(input_filename)
results <- as.data.frame(data$data)
colnames(results) <- map(data$labels, ~ .x[[1]])

# Log-scale the SNR to make the distribution closer to normal and use dB
results$LogSNR <- 10*log10(results$SNR)

results$Subject <- as.factor(results$Subject)
results$Measure <- as.factor(results$Measure)
results$Type <- as.factor(results$Type)
levels(results$Measure) <- unlist(data$measures)
levels(results$Type) <- unlist(data$types)


results_avg <- results %>% 
  aggregate(. ~ Subject + Measure + Type, mean)

results_snr_avg <- results_avg %>%
  aggregate(. ~ Subject, first)

results_snr <- results %>%
  aggregate(. ~ Subject + Session, first)


## Between-subject effect of SNR on accuracy
(with(results_snr_avg, cor.test(LogSNR, Accuracy)))

p_snr_acc_between <- ggplot(results_snr_avg, aes(x = LogSNR, y = Accuracy)) +
  geom_point() + 
  geom_smooth(method = 'lm', color = 'blue') +
  theme_classic()


## Within-subject effect of SNR on accuracy
fm_within <- scaleAndFitLME(results_snr, snr_vs_accuracy_formula, c('LogSNR', 'Accuracy'))
summary(fm_within)

# Restore original betas
beta0 <- lmerTest:::get_coefmat(fm_within)["(Intercept)","Estimate"]
beta1 <- lmerTest:::get_coefmat(fm_within)["LogSNR","Estimate"]
c(beta0_orig, beta1_orig) %<-% restoreBetas(beta0, beta1, results_snr, 
                                            x = 'LogSNR', y = 'Accuracy')

p_snr_acc_within <- ggplot(results_snr, aes(x = LogSNR, y = Accuracy)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(intercept = beta0_orig, slope = beta1_orig, 
              color = "blue", linewidth = 1) +
  scale_color_grey() +
  xlab('SNR (dB)') + ylab('Accuracy') +
  theme_classic()


## Longitudinal changes in SNR
fm_session <- scaleAndFitLME(results_snr, session_vs_snr_formula, c('Session', 'LogSNR'))
summary(fm_session)

# Restore original betas
beta0 <- lmerTest:::get_coefmat(fm_session)["(Intercept)","Estimate"]
beta1 <- lmerTest:::get_coefmat(fm_session)["Session","Estimate"]
c(beta0_orig, beta1_orig) %<-% restoreBetas(beta0, beta1, results_snr, 
                                            x = 'Session', y = 'LogSNR')

p_snr_session <- ggplot(results_snr, aes(x = Session, y = LogSNR)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(intercept = beta0_orig, slope = beta1_orig, 
              color = "blue", linewidth = 1) +
  scale_color_grey() +
  xlab('Session') + ylab('SNR (dB)') +
  theme_classic()


fig_snr <- plot_grid(p_snr_acc_between, p_snr_acc_within, p_snr_session, 
                     nrow = 1, labels = c('A', 'B', 'C'))
save_plot(file.path(output.folder, 'laplace-snr-rest-accuracy-session.png'), 
          fig_snr, base_width = 9, base_height = 3, bg = "white")


## SNR ~ Connectivity

laplace_snr_vs_connectivity_formula <- as.formula(
  str_replace_all(snr_vs_connectivity_formula_base, c("PS" = "Connectivity")))
laplace_conn_snr <- results %>%
  nest_by(Measure, Type) %>%
  mutate(fm = list(scaleAndFitLME(data, laplace_snr_vs_connectivity_formula, 
                                  c("LogSNR", "Connectivity")))) %>%
  mutate(intercept = list(lmerTest:::get_coefmat(fm)["(Intercept)",]),
         slope = list(lmerTest:::get_coefmat(fm)["LogSNR",]),
         ci = list(confint(fm))) %>%
  mutate(beta0 = intercept[["Estimate"]],
         beta1 = slope[["Estimate"]],
         ci_min = ci[["LogSNR", 1]],
         ci_max = ci[["LogSNR", 2]]) %>%
  mutate(betas_orig = list(restoreBetas(beta0, beta1, data, "LogSNR", "Connectivity"))) %>%
  mutate(beta0_orig = betas_orig[[1]], 
         beta1_orig = betas_orig[[2]])
laplace_conn_snr_df <- as.data.frame(do.call(rbind, 
                                             lapply(as.list(t(laplace_conn_snr$slope)), unlist)))
laplace_conn_snr_df$Significant <- laplace_conn_snr_df$`Pr(>|t|)` < p.threshold
laplace_conn_snr_df$Significant.MC <- laplace_conn_snr_df$`Pr(>|t|)` < p.mc
laplace_conn_snr_df$Predictor <- 'SNR'
laplace_conn_snr_df <- cbind(laplace_conn_snr[, c('Measure', 'Type', 'ci_min', 'ci_max')],
                             laplace_conn_snr_df)

dat_annot <- laplace_conn_snr_df
dat_annot$y <- c(0.31, 0.16, 0.38, 0.18, 0.95, 0.76)
p_snr_connectivity <- ggplot(results, aes(x = LogSNR, y = Connectivity)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(aes(intercept = beta0_orig, slope = beta1_orig),
              data = laplace_conn_snr, color = "blue") +
  geom_text(data = dat_annot,
            mapping = aes(x = -2, y = y, 
                          label = sapply(Estimate, \(x) paste0('β = ', format(x, digits = 1))),
                          fontface = sapply(`Pr(>|t|)`, \(x) mark(x))),
            size = 3.5, hjust = 0, vjust = 0) +
  facet_wrap(Type ~ Measure, scales = "free") +
  scale_color_grey() +
  xlab('SNR (dB)') + ylab('Connectivity') +
  theme_classic() +
  theme(strip.background = element_blank())

save_plot(file.path(output.folder, 'figCsupp1-laplace-snr-connectivity.png'), 
          p_snr_connectivity, base_width = 9, base_height = 6, bg = "white")


## Plot distributions

p_connectivity_dist <- ggplot(results, aes(x = Connectivity)) +
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 50, alpha = 0.5) +
  geom_density() +
  facet_wrap(Type ~ Measure, scales = "free") +
  theme_classic()

save_plot(file.path(output.folder, 'laplace-connectivity-density-plots.png'), 
          p_connectivity_dist, base_width = 9, base_height = 6, bg = "white")


## Between-subject effect of connectivity on performance

conn_laplace_avg <- results %>% 
  aggregate(. ~ Subject + Measure + Type, mean)


# Before correction for SNR
laplace_conn_acc_avg_corr <- conn_laplace_avg %>%
  group_by(Measure, Type) %>%
  do(c = cor.test(~ Connectivity + Accuracy, data = .))
laplace_conn_acc_avg_cor_data <- as.list(t(laplace_conn_acc_avg_corr[,'c']))
laplace_conn_acc_avg_corr_df <- as.data.frame(do.call(rbind, 
                                        lapply(laplace_conn_acc_avg_cor_data, unlist)))
laplace_conn_acc_avg_corr_df$Significant <- laplace_conn_acc_avg_corr_df$p.value < p.threshold
laplace_conn_acc_avg_corr_df$Significant.MC <- laplace_conn_acc_avg_corr_df$p.value < p.mc


p_conn_acc_between <- ggplot(conn_laplace_avg, aes(x = Connectivity, y = Accuracy)) +
  geom_point() +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_smooth(color = "blue", method = "lm") +
  facet_wrap(Type ~ Measure, scales = "free") +
  theme_classic()


save_plot(file.path(output.folder, 'laplace-connectivity-accuracy-between.png'), 
          p_conn_acc_between, base_width = 9, base_height = 6, bg = "white")

  
# After correction for SNR
laplace_conn_acc_avg_corr_SNR <- conn_laplace_avg %>%
  group_by(Measure, Type) %>%
  do(c = pcor.test(.data$Connectivity, .data$Accuracy, .data$LogSNR))
laplace_conn_acc_avg_corr_SNR_data <- as.list(t(laplace_conn_acc_avg_corr_SNR[,'c']))
laplace_conn_acc_avg_corr_SNR_df <- as.data.frame(do.call(rbind, 
                                                      laplace_conn_acc_avg_corr_SNR_data))
laplace_conn_acc_avg_corr_SNR_df$Significant <- laplace_conn_acc_avg_corr_SNR_df$p.value < p.threshold
laplace_conn_acc_avg_corr_SNR_df$Significant.MC <- laplace_conn_acc_avg_corr_SNR_df$p.value < p.mc


## Within-subject effect of connectivity on performance

laplace_conn_SNRcorrect <- results %>% 
  group_by(Subject, Measure, Type) %>%
  mutate(ConnSNR = lm(Connectivity ~ LogSNR)$residuals,
         AccSNR = lm(Accuracy ~ LogSNR)$residuals) %>%
  ungroup()

# Before correction for SNR
laplace_accuracy_vs_connectivity_formula <- as.formula(
  str_replace_all(accuracy_vs_connectivity_formula_base, c("PS" = "Connectivity")))
laplace_conn_acc_within <- results %>%
  nest_by(Measure, Type) %>%
  mutate(fm = list(scaleAndFitLME(data, laplace_accuracy_vs_connectivity_formula, 
                                  c("Connectivity", "Accuracy")))) %>%
  mutate(intercept = list(lmerTest:::get_coefmat(fm)["(Intercept)",]),
         slope = list(lmerTest:::get_coefmat(fm)["Connectivity",]),
         ci = list(confint(fm))) %>%
  mutate(beta0 = intercept[["Estimate"]],
         beta1 = slope[["Estimate"]],
         ci_min = ci[["Connectivity", 1]],
         ci_max = ci[["Connectivity", 2]]) %>%
  mutate(betas_orig = list(restoreBetas(beta0, beta1, data, "Connectivity", "Accuracy"))) %>%
  mutate(beta0_orig = betas_orig[[1]], 
         beta1_orig = betas_orig[[2]])
laplace_conn_acc_within_df <- as.data.frame(do.call(rbind, 
                                             lapply(as.list(t(laplace_conn_acc_within$slope)), unlist)))
laplace_conn_acc_within_df$Significant <- laplace_conn_acc_within_df$`Pr(>|t|)` < p.threshold
laplace_conn_acc_within_df$Significant.MC <- laplace_conn_acc_within_df$`Pr(>|t|)` < p.mc
laplace_conn_acc_within_df$Response <- 'Accuracy'
laplace_conn_acc_within_df <- cbind(laplace_conn_acc_within[, c('Measure', 'Type', 'ci_min', 'ci_max')],
                                    laplace_conn_acc_within_df)

dat_annot <- laplace_conn_acc_within_df
dat_annot$x <- c(0.32, 0.17, 0.38, 0.18, 0.95, 0.77)
p_conn_acc_within <- ggplot(results, aes(x = Connectivity, y = Accuracy)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(aes(intercept = beta0_orig, slope = beta1_orig),
              data = laplace_conn_acc_within, color = "blue") +
  geom_text(data = dat_annot,
            mapping = aes(x = x, y = 0.25, 
                          label = sapply(Estimate, \(x) paste0('β = ', format(x, digits = 1))),
                          fontface = sapply(`Pr(>|t|)`, \(x) mark(x))),
            size = 3.5, hjust = 1, vjust = 0) +
  facet_wrap(Type ~ Measure, scales = "free") +
  scale_color_grey() +
  xlab('Connectivity') + ylab('Accuracy') +
  theme_classic() +
  theme(strip.background = element_blank())


save_plot(file.path(output.folder, 'figCsupp2-laplace-connectivity-performance-within.png'), 
          p_conn_acc_within, base_width = 9, base_height = 6, bg = "white")


# After correction for SNR
laplace_accuracy_vs_connectivity_formula_SNRcorr <- as.formula(
  str_replace_all(accuracy_vs_connectivity_formula_SNRCorr, c("PS" = "Connectivity")))
laplace_conn_acc_within_SNR <- results %>%
  nest_by(Measure, Type) %>%
  mutate(fm = list(scaleAndFitLME(data, laplace_accuracy_vs_connectivity_formula_SNRcorr, 
                                  c("Connectivity", "Accuracy")))) %>%
  mutate(intercept = list(lmerTest:::get_coefmat(fm)["(Intercept)",]),
         slope = list(lmerTest:::get_coefmat(fm)["Connectivity",]),
         ci = list(confint(fm))) %>%
  mutate(beta0 = intercept[["Estimate"]],
         beta1 = slope[["Estimate"]],
         ci_min = ci[["Connectivity", 1]],
         ci_max = ci[["Connectivity", 2]]) %>%
  mutate(betas_orig = list(restoreBetas(beta0, beta1, data, "Connectivity", "Accuracy"))) %>%
  mutate(beta0_orig = betas_orig[[1]], 
         beta1_orig = betas_orig[[2]])
laplace_conn_acc_within_SNR_df <- as.data.frame(do.call(rbind, 
                                                    lapply(as.list(t(laplace_conn_acc_within_SNR$slope)), unlist)))
laplace_conn_acc_within_SNR_df$Significant <- laplace_conn_acc_within_SNR_df$`Pr(>|t|)` < p.threshold
laplace_conn_acc_within_SNR_df$Significant.MC <- laplace_conn_acc_within_SNR_df$`Pr(>|t|)` < p.mc
laplace_conn_acc_within_SNR_df$Response <- 'Accuracy $|$ SNR'
laplace_conn_acc_within_SNR_df <- cbind(laplace_conn_acc_within_SNR[, c('Measure', 'Type', 'ci_min', 'ci_max')],
                                        laplace_conn_acc_within_SNR_df)


## Longitudinal changes in connectivity

# Before correction for SNR
laplace_connectivity_vs_session_formula <- as.formula(
  str_replace_all(connectivity_vs_session_formula_base, c("PS" = "Connectivity")))
laplace_conn_session_within <- results %>%
  nest_by(Measure, Type) %>%
  mutate(fm = list(scaleAndFitLME(data, laplace_connectivity_vs_session_formula, 
                                  c("Session", "Connectivity")))) %>%
  mutate(intercept = list(lmerTest:::get_coefmat(fm)["(Intercept)",]),
         slope = list(lmerTest:::get_coefmat(fm)["Session",]),
         ci = list(confint(fm))) %>%
  mutate(beta0 = intercept[["Estimate"]],
         beta1 = slope[["Estimate"]],
         ci_min = ci[["Session", 1]],
         ci_max = ci[["Session", 2]]) %>%
  mutate(betas_orig = list(restoreBetas(beta0, beta1, data, "Session", "Connectivity"))) %>%
  mutate(beta0_orig = betas_orig[[1]], 
         beta1_orig = betas_orig[[2]])
laplace_conn_session_within_df <- as.data.frame(do.call(rbind, 
                                                    lapply(as.list(t(laplace_conn_session_within$slope)), unlist)))
laplace_conn_session_within_df$Significant <- laplace_conn_session_within_df$`Pr(>|t|)` < p.threshold
laplace_conn_session_within_df$Significant.MC <- laplace_conn_session_within_df$`Pr(>|t|)` < p.mc
laplace_conn_session_within_df$Predictor <- 'Session'
laplace_conn_session_within_df <- cbind(laplace_conn_session_within[, c('Measure', 'Type', 'ci_min', 'ci_max')],
                                        laplace_conn_session_within_df)


dat_annot <- laplace_conn_session_within_df
dat_annot$y <- c(0.33, 0.16, 0.4, 0.19, 0.95, 0.79)
p_conn_session_within <- ggplot(results, aes(x = Session, y = Connectivity)) +
  geom_point(mapping = aes(color = Subject), size = 1, show.legend = F) +
  geom_smooth(mapping = aes(color = Subject, group = Subject), linetype = 1, 
              linewidth = 0.5, method = 'lm', se = F, show.legend = F) +
  geom_abline(aes(intercept = beta0_orig, slope = beta1_orig),
              data = laplace_conn_session_within, color = "blue") +
  geom_text(data = dat_annot,
            mapping = aes(x = 11, y = y, 
                          label = sapply(Estimate, \(x) paste0('β = ', format(x, digits = 1))),
                          fontface = sapply(`Pr(>|t|)`, \(x) mark(x))),
            size = 3.5, hjust = 1, vjust = 1) +
  facet_wrap(Type ~ Measure, scales = "free") +
  scale_color_grey() +
  xlab('Session') + ylab('Connectivity') +
  theme_classic() +
  theme(strip.background = element_blank())


save_plot(file.path(output.folder, 'figCsupp3-laplace-connectivity-longitude.png'), 
          p_conn_session_within, base_width = 9, base_height = 6, bg = "white")


# After correction for SNR
laplace_connectivity_vs_session_formula_SNRcorr <- as.formula(
  str_replace_all(connectivity_vs_session_formula_SNRCorr, c("PS" = "Connectivity")))
laplace_conn_session_within_SNR <- results %>%
  nest_by(Measure, Type) %>%
  mutate(fm = list(scaleAndFitLME(data, laplace_connectivity_vs_session_formula_SNRcorr, 
                                  c("Session", "Connectivity")))) %>%
  mutate(intercept = list(lmerTest:::get_coefmat(fm)["(Intercept)",]),
         slope = list(lmerTest:::get_coefmat(fm)["Session",]),
         ci = list(confint(fm))) %>%
  mutate(beta0 = intercept[["Estimate"]],
         beta1 = slope[["Estimate"]],
         ci_min = ci[["Session", 1]],
         ci_max = ci[["Session", 2]]) %>%
  mutate(betas_orig = list(restoreBetas(beta0, beta1, data, "Session", "Connectivity"))) %>%
  mutate(beta0_orig = betas_orig[[1]], 
         beta1_orig = betas_orig[[2]])
laplace_conn_session_within_SNR_df <- as.data.frame(do.call(rbind, 
                                                        lapply(as.list(t(laplace_conn_session_within_SNR$slope)), unlist)))
laplace_conn_session_within_SNR_df$Significant <- laplace_conn_session_within_SNR_df$`Pr(>|t|)` < p.threshold
laplace_conn_session_within_SNR_df$Significant.MC <- laplace_conn_session_within_SNR_df$`Pr(>|t|)` < p.mc
laplace_conn_session_within_SNR_df$Predictor <- 'Session'
laplace_conn_session_within_SNR_df <- cbind(laplace_conn_session_within_SNR[, c('Measure', 'Type', 'ci_min', 'ci_max')],
                                            laplace_conn_session_within_SNR_df)

# Export the results to TeX
all_results <- bind_rows(
  laplace_conn_snr_df,
  laplace_conn_acc_within_df,
  laplace_conn_acc_within_SNR_df,
  laplace_conn_session_within_df,
  laplace_conn_session_within_SNR_df
) %>%
  mutate_at(vars('Type'), ~ case_when(
    . == 'Within-hemisphere' ~ 'WH',
    . == 'Across-hemisphere' ~ 'AH')
    ) %>%
  mutate(TypeMeasure = paste(Type, Measure))
names(all_results) <- c('Measure', 'Type', 'CIMin', 'CIMax', '$\\beta$', 'SE', 
                        'df', 't-value', 'p-value', 'Significant', 
                        'Significant.MC', 'Predictor', 'Response', 'TypeMeasure')
all_results$CIDisp <- paste(
  "$[", format(round(all_results$CIMin, 3), nsmall = 3), ",", 
  format(round(all_results$CIMax, 3), nsmall = 3), "]$", sep = "")
all_results$Response[1:6] <- all_results$TypeMeasure[1:6]
all_results$Predictor[7:18] <- all_results$TypeMeasure[7:18]
all_results$Response[19:24] <- all_results$TypeMeasure[19:24]
all_results$Response[25:30] <- paste(all_results$TypeMeasure[25:30], '$|$ SNR')

# Reorder rows to match joint multiverse analysis
new_order <- c(1, 3, 5, 2, 4, 6, 
               7, 9, 11, 13, 15, 17, 8, 10, 12, 14, 16, 18, 
               19, 21, 23, 25, 27, 29, 20, 22, 24, 26, 28, 30)

all_results_latex <- all_results[new_order, c('Predictor', 'Response',  
                                              '$\\beta$', 't-value', 
                                              'p-value', 'CIDisp')] %>%
  rename("95\\% CI" = "CIDisp") %>%
  mutate_at(vars('t-value', '$\\beta$', 'p-value'), round, 3)
sig.mc <- all_results[new_order,] %>%
  mutate_at(vars(Significant.MC), ~ if_else(., '*', '', missing = ''))

all_results_latex$`p-value` <- cell_spec(
  paste(if_else(all_results_latex[,'p-value'] < 0.001, '$<$0.001', 
                format(all_results_latex[,'p-value'], digits = 3)), 
        sig.mc[,'Significant.MC'], sep = ''), 
  format = 'latex', bold = all_results[new_order, 'Significant'], escape = F)

all_results_latex %>%
  knitr::kable(format = "latex", booktabs = T, escape = F, align = 'c', row.names = F,
               linesep = c(
                 "", "", "\\addlinespace", "", "", "\\midrule",
                 "", "", "\\addlinespace", "", "", "\\addlinespace", 
                 "", "", "\\addlinespace", "", "", "\\midrule",
                 "", "", "\\addlinespace", "", "", "\\addlinespace", 
                 "", "", "\\addlinespace", "", ""
               )) %>%
  tex.save(output_filename, "Summary", ., prefix = 'laplaceConnectivity')


### Save all the results
save(laplace_conn_snr_df, laplace_conn_acc_within_df, laplace_conn_acc_within_SNR_df,
     laplace_conn_session_within_df, laplace_conn_session_within_SNR_df,
     file = file.path(r.path, "laplace_connectivity.RData"))