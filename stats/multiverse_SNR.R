###
# Effects of SNR on accuracy and connectivity as well as changes over time in 
# the multiverse analysis
#
# Outputs:
# 1. fig4-snr-rest-accuracy-session.png
# 2. fig5-snr-connectivity.png
# 3. multiverse_SNR.RData - intermediate results of the analysis
###

# Load the data
c(results, results_SNR, pipeline_desc) %<-% load_multiverse(input_filename)

# Create output folder for images
multiverse.folder <- paste('../results/stats/', prefix, sep='')
dir.create(multiverse.folder)


### Between-subject effect of SNR on accuracy

# Split analysis - fit mixed models
SNR_vs_acc_avg_stats <- fitMultiverseSplitBetween(
  df = results_SNR,
  list(x = 'LogSNR', y = 'Accuracy', group = 'Subject',
       measure.name = '', measure.type = ''),
  pipeline_desc = pipeline_desc
) %>% dplyr::select(-Measure, -Type)

# Joint results - use data from all pipelines together
data_avg <- results_SNR %>%
  group_by(Pipeline, Subject) %>%
  summarise(LogSNR_avg = mean(LogSNR), Accuracy_avg = mean(Accuracy)) %>%
  ungroup()
SNR_vs_acc_between_joint <- rmcorr(Pipeline, LogSNR_avg, Accuracy_avg, 
                                   dataset = data_avg)
print(SNR_vs_acc_between_joint)

SNR_vs_accuracy_between_results <- list(
  predictor = 'Average LogSNR', response = 'Average Accuracy', level = 'Between',
  split_significant = sum(SNR_vs_acc_avg_stats$Significant),
  split_total = length(SNR_vs_acc_avg_stats$Significant),
  split_t.value = mean(SNR_vs_acc_avg_stats$t.value),
  joint_p.value = SNR_vs_acc_between_joint$p,
  joint_result = SNR_vs_acc_between_joint$p <= p.threshold,
  joint_result.MC = SNR_vs_acc_between_joint$p <= p.threshold  # No correction
)


### Plot the relationship between SNR and accuracy for pipelines separately
p_separate <- ggplot(data = data_avg) +
  geom_point(mapping = aes(x = LogSNR_avg, y = Accuracy_avg), size = 0.5) +
  geom_line(mapping = aes(x = LogSNR_avg, y = Accuracy_avg), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 2) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")
ggsave(file.path(plot.path, prefix,
                 paste(prefix, '_snr_vs_accuracy_pipelines_between.png', sep = '')),
       plot = p_separate, width = 7, height = 5)


### Within-subject Effect of SNR on Accuracy

# Split Analysis - fit mixed models
SNR_vs_acc_stats <- fitMultiverseSplitLME(
  list(formula_split = snr_vs_accuracy_formula,
       cols.scale = c('LogSNR', 'Accuracy'), predictor = 'LogSNR',
       measure.name = '', measure.type = ''), 
  df = results_SNR, pipeline_desc = pipeline_desc
) %>% dplyr::select(-Measure, -Type)
assert("Accuracy ~ ROI SNR [split] did not converge",
       all(SNR_vs_acc_stats$Converged))

# Joint Analysis - use the data from all pipelines together
SNR_vs_acc_joint <- scaleAndFitLME(
  results_SNR, formula = snr_vs_accuracy_formula_joint,
  columns_to_scale = c('LogSNR', 'Accuracy'))
assert("Accuracy ~ ROI SNR [joint] did not converge",
       has_converged(SNR_vs_acc_joint))
SNR_vs_acc_joint_ci <- confint(SNR_vs_acc_joint)
summary(SNR_vs_acc_joint)
print(SNR_vs_acc_joint_ci)
SNR_vs_acc_joint_stats <- lmerTest:::get_coefmat(SNR_vs_acc_joint)['LogSNR',] %>%
  t() %>% as.data.frame() %>%
  rename("p.value" = "Pr(>|t|)", "t.value" = "t value")

SNR_vs_accuracy_within_results <- list(
  predictor = 'SNR', response = 'Accuracy', level = 'Within', 
  split_significant = sum(SNR_vs_acc_stats$Significant),
  split_total = length(SNR_vs_acc_stats$Significant),
  split_t.value = mean(SNR_vs_acc_stats$t.value),
  joint_estimate = SNR_vs_acc_joint_stats$Estimate,
  joint_df = SNR_vs_acc_joint_stats$df,
  joint_t.value = SNR_vs_acc_joint_stats$t.value,
  joint_p.value = SNR_vs_acc_joint_stats$p.value,
  joint_ci_min = SNR_vs_acc_joint_ci[["LogSNR", 1]],
  joint_ci_max = SNR_vs_acc_joint_ci[["LogSNR", 2]],
  joint_result = SNR_vs_acc_joint_stats$p.value <= p.threshold,
  joint_result.MC = SNR_vs_acc_joint_stats$p.value <= p.threshold  # No correction
)


### Plot within-subject effect for different pipelines
p_separate <- ggplot(data = results_SNR) +
  geom_point(mapping = aes(x = LogSNR, y = Accuracy, group = Subject, color = Subject), size = 0.5) +
  geom_line(mapping = aes(x = LogSNR, y = Accuracy, group = Subject), color = "black", alpha = 0.5, 
            stat="smooth", method = "lm", se = F, linewidth = 0.5) +
  geom_line(mapping = aes(x = LogSNR, y = Accuracy), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 2) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")
ggsave(file.path(plot.path, prefix,
                 paste(prefix, '_snr_vs_accuracy_pipelines_within.png', sep = '')),
       plot = p_separate, width = 7, height = 5)


### Longitudinal Changes in SNR

# Split Analysis - fit mixed models
SNR_vs_session_stats <- fitMultiverseSplitLME(
  list(formula_split = session_vs_snr_formula,
       cols.scale = c('LogSNR', 'Session'), predictor = 'Session',
       measure.name = '', measure.type = ''), 
  df = results_SNR, pipeline_desc = pipeline_desc
) %>% dplyr::select(-Measure, -Type)
assert("ROI SNR ~ Session [split] did not converge",
       all(SNR_vs_session_stats$Converged))

# Joint Analysis - use the data from all pipelines together
SNR_vs_session_joint <- scaleAndFitLME(
  results_SNR, formula = session_vs_snr_formula_joint,
  columns_to_scale = c('LogSNR', 'Session'))
assert("ROI SNR ~ Session [joint] did not converge",
       has_converged(SNR_vs_session_joint))
SNR_vs_session_joint_ci <- confint(SNR_vs_session_joint)
summary(SNR_vs_session_joint)
print(SNR_vs_session_joint_ci)
SNR_vs_session_joint_stats <- lmerTest:::get_coefmat(SNR_vs_session_joint)['Session',] %>%
  t() %>% as.data.frame() %>%
  rename("p.value" = "Pr(>|t|)", "t.value" = "t value")

SNR_vs_session_results <- list(
  predictor = 'Session', response = 'SNR', level = 'Within', 
  split_significant = sum(SNR_vs_session_stats$Significant),
  split_total = length(SNR_vs_session_stats$Significant),
  split_t.value = mean(SNR_vs_session_stats$t.value),
  joint_estimate = SNR_vs_session_joint_stats$Estimate,
  joint_df = SNR_vs_session_joint_stats$df,
  joint_t.value = SNR_vs_session_joint_stats$t.value,
  joint_p.value = SNR_vs_session_joint_stats$p.value,
  joint_ci_min = SNR_vs_session_joint_ci[["Session", 1]],
  joint_ci_max = SNR_vs_session_joint_ci[["Session", 2]],
  joint_result = SNR_vs_session_joint_stats$p.value <= p.threshold,
  joint_result.MC = SNR_vs_session_joint_stats$p.value <= p.threshold  # No correction
)
print(SNR_vs_session_results)


### Plot the multiverse for SNR
lim_actual <- max(abs(c(SNR_vs_acc_avg_stats$Estimate, 
                        SNR_vs_acc_stats$Estimate, 
                        SNR_vs_session_stats$Estimate)))
lim_snr_acc = ceilingn(lim_actual, 1)

# BB -> Broadband
levels(SNR_vs_acc_avg_stats$Band) <- c('Broadband', 'Narrowband')
levels(SNR_vs_acc_stats$Band) <- c('Broadband', 'Narrowband')
levels(SNR_vs_session_stats$Band) <- c('Broadband', 'Narrowband')

p1 <- plotMultiverseSplit(SNR_vs_acc_avg_stats, 'Estimate', 'Significant',
                          lim = lim_snr_acc, facet_rule = 'Band ~ Mask')

p2 <- plotMultiverseSplit(SNR_vs_acc_stats, 'Estimate', 'Significant',
                          lim = lim_snr_acc, facet_rule = 'Band ~ Mask')

p3 <- plotMultiverseSplit(SNR_vs_session_stats, 'Estimate', 'Significant',
                          lim = lim_snr_acc, facet_rule = 'Band ~ Mask')

p_multiverse = list(between = p1, within = p2, session = p3)


### Add labels to distinguish between Laplace and ROI SNR
p_label_laplace = ggdraw() +
  draw_label("Average SNR of the mu rhythm measured from C3- and C4-Laplace",
             fontface = 'bold', y = 0.75)
p_label_roi = ggdraw() +
  draw_label("Average SNR of the mu rhythm in the sensorimotor ROIs",
             fontface = 'bold', y = 0.25, vjust = 0)

### Plot the effects of Laplace and ROI SNR on Accuracy -- Figure 4
fig4_step1 <- plot_grid(
  p_laplace$examples, p_laplace$group_cmp, p_laplace$dynamics,
  p_laplace$between, p_laplace$within, p_laplace$session,
  nrow = 2, ncol = 3, align = 'hv', axis = 'l',
  rel_heights = c(1.5, 2),
  labels = 'AUTO'
)

fig4_step2 <- plot_grid(
  p_multiverse$between + theme(legend.position = 'none'),
  p_multiverse$within + theme(legend.position = 'none'),
  p_multiverse$session + theme(legend.position = 'none'),
  ncol = 3, labels = c('G', 'H', 'I')
)

fig4_step3 <- plot_grid(
  p_label_laplace, 
  fig4_step1, 
  p_label_roi, 
  fig4_step2, 
  get_legend(p_multiverse$between),
  ncol = 1, rel_heights = c(0.25, 3.5, 0.4, 1.4, 0.82),
  labels = c('', '', '', '')
)

save_plot(file.path(plot.path, prefix, 'fig4-snr-rest-accuracy-session.pdf'), 
          fig4_step3, base_width = 9, base_height = 8)
save_plot(file.path(plot.path, prefix, 'fig4-snr-rest-accuracy-session.png'), 
          fig4_step3, bg = "white", base_width = 9, base_height = 8)


### Effect of SNR on Connectivity
param_space <- list(
  list(formula_split = snr_vs_connectivity_formula("ImCoh_Within"),
       formula_joint = snr_vs_connectivity_formula("ImCoh_Within", joint = T),
       cols.scale = c('LogSNR', 'ImCoh_Within'), 
       predictor = 'LogSNR', response = 'WH ImCoh',
       measure.name = 'ImCoh', measure.type = 'Within'),
  list(formula_split = snr_vs_connectivity_formula("LagCoh_Within"),
       formula_joint = snr_vs_connectivity_formula("LagCoh_Within", joint = T),
       cols.scale = c('LogSNR', 'LagCoh_Within'), 
       predictor = 'LogSNR', response = 'WH LagCoh',
       measure.name = 'LagCoh', measure.type = 'Within'),
  list(formula_split = snr_vs_connectivity_formula("Coh_Within"),
       formula_joint = snr_vs_connectivity_formula("Coh_Within", joint = T),
       cols.scale = c('LogSNR', 'Coh_Within'), 
       predictor = 'LogSNR', response = 'WH Coherence',
       measure.name = 'Coherence', measure.type = 'Within'),
  list(formula_split = snr_vs_connectivity_formula("ImCoh_Across"),
       formula_joint = snr_vs_connectivity_formula("ImCoh_Across", joint = T),
       cols.scale = c('LogSNR', 'ImCoh_Across'), 
       predictor = 'LogSNR', response = 'AH ImCoh',
       measure.name = 'ImCoh', measure.type = 'Across'),
  list(formula_split = snr_vs_connectivity_formula("LagCoh_Across"),
       formula_joint = snr_vs_connectivity_formula("LagCoh_Across", joint = T),
       cols.scale = c('LogSNR', 'LagCoh_Across'), 
       predictor = 'LogSNR', response = 'AH LagCoh',
       measure.name = 'LagCoh', measure.type = 'Across'),
  list(formula_split = snr_vs_connectivity_formula("Coh_Across"),
       formula_joint = snr_vs_connectivity_formula("Coh_Across", joint = T),
       cols.scale = c('LogSNR', 'Coh_Across'),
       predictor = 'LogSNR', response = 'AH Coherence',
       measure.name = 'Coherence', measure.type = 'Across')
)

snr_vs_conn_stats <- lapply(param_space, fitMultiverseSplitLME, 
                            df = results, pipeline_desc = pipeline_desc)
snr_vs_conn_stats <- do.call(rbind, snr_vs_conn_stats)
snr_vs_conn_stats$Significant.MC <- snr_vs_conn_stats$p.value < p.mc
assert("PS ~ ROI SNR [split] did not converge",
       all(snr_vs_conn_stats$Converged))

snr_vs_conn_joint_stats <- lapply(param_space, fitMultiverseJointLME, 
                                  df = results)
snr_vs_conn_joint_stats <- do.call(rbind, snr_vs_conn_joint_stats)
assert("PS ~ ROI SNR [joint] did not converge",
       all(snr_vs_conn_joint_stats$Converged))

SNR_vs_conn_results <- lapply(param_space, getSummary,
                              split_stats = snr_vs_conn_stats,
                              joint_stats = snr_vs_conn_joint_stats,
                              commonFields = list(level = 'Within'))


### Multiverse Plot for SNR vs Connectivity --- Figure 5
snr_vs_conn_stats$fMeasure = factor(snr_vs_conn_stats$Measure, 
                                    levels = c('ImCoh', 'LagCoh', 'Coherence'))
snr_vs_conn_stats$fType = factor(snr_vs_conn_stats$Type, 
                                 levels = c('Within', 'Across'))
levels(snr_vs_conn_stats$fType) <- list("Within-hemisphere" = "Within",
                                        "Across-hemisphere" = "Across")
levels(snr_vs_conn_stats$Band) <- list("Broadband" = "BB", "Narrowband" = "NB")

p_snr_conn <- plotMultiverseSplit(snr_vs_conn_stats, 'Estimate', 'Significant.MC',
                                  facet_rule = 'fType + Band ~ fMeasure + Mask')


### Sanity check: spectra of different connectivity measures
fdispmin <- 3
fdispmax <- 40
colors <- c("eLORETA" = "#1984c5")
linetypes <- c("1SVD" = "solid", "3SVD" = "42", "AVG-F" = "11")
selection <- (conn_df$Inverse == "eLORETA" & conn_df$Mask == "No Mask" & 
  conn_df$Freqs >= fdispmin & conn_df$Freqs <= fdispmax)
p_imcoh <- ggplot(conn_df[selection & conn_df$Measure == "ImCoh",]) + 
  geom_line(aes(x = Freqs, y = Conn, color = Inverse, 
                linetype = ROI_Method, group = Pipeline), 
            linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  xlim(fdispmin, fdispmax) +
  xlab('') + ylab('Phase Synchronization') +
  facet_nested(Type ~ Measure, scales = 'free', switch = 'y') +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(face = 'bold'))

p_lagcoh <- ggplot(conn_df[selection & conn_df$Measure == "LagCoh",]) + 
  geom_line(aes(x = Freqs, y = Conn, color = Inverse,
                linetype = ROI_Method, group = Pipeline),
            linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  xlim(fdispmin, fdispmax) +
  xlab('Frequency (Hz)') +
  facet_grid(Type ~ Measure, scales = 'free', switch = 'y') +
  theme_classic() +
  theme(axis.title.y.left = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'bold'),
        strip.text.y = element_blank())

p_coh <- ggplot(conn_df[selection & conn_df$Measure == "Coherence",]) + 
  geom_line(aes(x = Freqs, y = Conn, color = Inverse,
                linetype = ROI_Method, group = Pipeline), 
            linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  xlim(fdispmin, fdispmax) +
  xlab('') +
  facet_grid(Type ~ Measure, scales = 'free', switch = 'y') +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'bold'),
        strip.text.y = element_blank())

### Combine the plots
legend <- get_legend(p_imcoh +
  guides(linetype = guide_legend(title = "ROI Aggregation Method"),
         color = guide_legend(title = 'Inverse Method')) +
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom',
        legend.key.width = unit(1, "cm")))
fig5A <- plot_grid(p_imcoh + theme(legend.position = 'none'),
                   p_lagcoh + theme(legend.position = 'none'),
                   p_coh + theme(legend.position = 'none'),
                   nrow = 1, rel_widths = c(0.35, 0.3, 0.35),
                   labels = c('A', '', ''))
fig5 <- plot_grid(fig5A, legend, p_snr_conn, 
                  ncol = 1, rel_heights = c(0.4, 0.05, 0.6),
                  labels = c('', '', 'B'))
save_plot(file.path(plot.path, prefix, 'fig5-snr-connectivity.pdf'),
          plot = fig5, base_width = 9, base_height = 11)
save_plot(file.path(plot.path, prefix, 'fig5-snr-connectivity.png'),
          plot = fig5, bg = "white", base_width = 9, base_height = 11)


### Plot SNR vs Connectivity for different pipelines

# ImCoh_Within
p_imcoh_within <- ggplot(data = results) +
  geom_point(mapping = aes(x = LogSNR, y = ImCoh_Within, group = Subject, color = Subject), size = 0.5) +
  geom_line(mapping = aes(x = LogSNR, y = ImCoh_Within, group = Subject), color = "black", alpha = 0.5, 
            stat="smooth", method = "lm", se = F, linewidth = 0.5) +
  geom_line(mapping = aes(x = LogSNR, y = ImCoh_Within), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 4) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")
ggsave(file.path(plot.path, prefix,
                 paste(prefix, '_snr_imcoh_within_pipelines.png', sep = '')),
       plot = p_imcoh_within, width = 7, height = 5)


# Coh_Within
p_coh_within <- ggplot(data = results) +
  geom_point(mapping = aes(x = LogSNR, y = Coh_Within, group = Subject, color = Subject), size = 0.5) +
  geom_line(mapping = aes(x = LogSNR, y = Coh_Within, group = Subject), color = "black", alpha = 0.5, 
            stat="smooth", method = "lm", se = F, linewidth = 0.5) +
  geom_line(mapping = aes(x = LogSNR, y = Coh_Within), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 4) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")
ggsave(file.path(plot.path, prefix,
                 paste(prefix, '_snr_coh_within_pipelines.png', sep = '')),
       plot = p_coh_within, width = 7, height = 5)


### Save all the results
save(SNR_vs_acc_avg_stats, SNR_vs_acc_between_joint, 
     SNR_vs_acc_stats, SNR_vs_acc_joint, SNR_vs_acc_joint_ci, SNR_vs_acc_joint_stats, 
     SNR_vs_session_stats, SNR_vs_session_joint, SNR_vs_session_joint_ci, SNR_vs_session_joint_stats,
     snr_vs_conn_stats, snr_vs_conn_joint_stats,
     SNR_vs_accuracy_between_results, SNR_vs_accuracy_within_results, 
     SNR_vs_session_results, SNR_vs_conn_results,
     results_SNR, data_avg,
     file = file.path(r.path, 'multiverse_SNR.RData'))
