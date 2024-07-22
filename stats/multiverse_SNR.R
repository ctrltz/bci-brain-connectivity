###
# Effects of SNR on accuracy and longitudinal changes in SNR in the multiverse
#
# Outputs:
# 1. fig5-snr-rest-accuracy-session.png
# 2. multiverse_SNR.RData - intermediate results of the analysis
###

# Load the data
c(results, pipeline_desc) %<-% load_multiverse(input_filename)

# Create output folder for images
prefix = "multiverse_SNR"
multiverse.folder <- paste('../results/stats/', prefix, sep='')
dir.create(multiverse.folder)


### Check for outliers
assert("Checking for ROI SNR outliers [split]",
       all(rbind(by(results, results$Pipeline, 
                    \(x) sum(checkOutliers(x$LogSNR)))) == 0))
assert("Checking for ROI SNR outliers [joint]",
       sum(checkOutliers(results$LogSNR)) == 0)


### Between-subject effect of SNR on accuracy

# Split analysis - fit mixed models
SNR_vs_acc_avg_stats <- fitMultiverseSplitBetween(
  df = results,
  list(x = 'LogSNR', y = 'Accuracy', group = 'Subject',
       measure.name = '', measure.type = ''),
  pipeline_desc = pipeline_desc
) %>% dplyr::select(-Measure, -Type)

# Joint results - use data from all pipelines together
data_avg <- results %>%
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
p_snr_acc_separate <- plotMultiverseScatterBetween(data_avg, 'LogSNR_avg', 'Accuracy_avg')
ggsave(file.path(plot.path, prefix,
                 paste('snr_vs_accuracy_pipelines_between.png', sep = '')),
       plot = p_snr_acc_separate, width = 9, height = 6)


### Within-subject Effect of SNR on Accuracy

# Split Analysis - fit mixed models
SNR_vs_acc_stats <- fitMultiverseSplitLME(
  list(formula_split = snr_vs_accuracy_formula,
       cols.scale = c('LogSNR', 'Accuracy'), predictor = 'LogSNR',
       measure.name = '', measure.type = ''), 
  df = results, pipeline_desc = pipeline_desc
) %>% dplyr::select(-Measure, -Type)
assert("Accuracy ~ ROI SNR [split] did not converge",
       all(SNR_vs_acc_stats$Converged))

# Joint Analysis - use the data from all pipelines together
SNR_vs_acc_joint <- scaleAndFitLME(
  results, formula = snr_vs_accuracy_formula_joint,
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
  split_significant_samedir = sum(SNR_vs_acc_stats$Significant &
                                  (SNR_vs_acc_stats$t.value * SNR_vs_acc_joint_stats$t.value > 0)),
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
p_snr_acc_separate <- plotMultiverseScatterWithin(results, 'LogSNR', 'Accuracy')
ggsave(file.path(plot.path, prefix,
                 paste('snr_vs_accuracy_pipelines_within.png', sep = '')),
       plot = p_snr_acc_separate, width = 9, height = 6)


### Longitudinal Changes in SNR

# Split Analysis - fit mixed models
SNR_vs_session_stats <- fitMultiverseSplitLME(
  list(formula_split = session_vs_snr_formula,
       cols.scale = c('LogSNR', 'Session'), predictor = 'Session',
       measure.name = '', measure.type = ''), 
  df = results, pipeline_desc = pipeline_desc, refit = T
) %>% dplyr::select(-Measure, -Type)
assert("ROI SNR ~ Session [split] did not converge",
       all(SNR_vs_session_stats$Converged))

# Joint Analysis - use the data from all pipelines together
SNR_vs_session_joint <- scaleAndFitLME(
  results, formula = session_vs_snr_formula_joint,
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
  split_significant_samedir = sum(SNR_vs_session_stats$Significant &
                                  (SNR_vs_session_stats$t.value * SNR_vs_session_joint_stats$t.value > 0)),
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


### Plot the longitudinal changes in SNR for pipelines separately
p_snr_session_separate <- plotMultiverseScatterWithin(results, 'Session', 'LogSNR')
ggsave(file.path(plot.path, prefix,
                 paste('snr_vs_session_pipelines_within.png', sep = '')),
       plot = p_snr_session_separate, width = 9, height = 6)


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
  draw_label("Average SNR of the mu rhythm in the sensorimotor ROIs (source space)",
             fontface = 'bold', y = 0.25, vjust = 0)

### Plot the effects of Laplace and ROI SNR on Accuracy -- Figure 5
fig5_step1 <- plot_grid(
  p_laplace$examples, p_laplace$group_cmp, p_laplace$dynamics,
  p_laplace$between, p_laplace$within, p_laplace$session,
  nrow = 2, ncol = 3, align = 'hv', axis = 'l',
  rel_heights = c(1.5, 2),
  labels = 'AUTO'
)

fig5_step2 <- plot_grid(
  p_multiverse$between + theme(legend.position = 'none'),
  p_multiverse$within + theme(legend.position = 'none'),
  p_multiverse$session + theme(legend.position = 'none'),
  ncol = 3, labels = c('G', 'H', 'I')
)

fig5_step3 <- plot_grid(
  p_label_laplace, 
  fig5_step1, 
  p_label_roi, 
  fig5_step2, 
  get_legend(p_multiverse$between),
  ncol = 1, rel_heights = c(0.25, 3.5, 0.4, 2.2, 0.82),
  labels = c('', '', '', '')
)

save_plot(file.path(plot.path, prefix, 'fig5-snr-rest-accuracy-session.pdf'), 
          fig5_step3, base_width = 9, base_height = 9.5)
save_plot(file.path(plot.path, prefix, 'fig5-snr-rest-accuracy-session.png'), 
          fig5_step3, bg = "white", base_width = 9, base_height = 9.5)


### Save all the results
save(SNR_vs_acc_avg_stats, SNR_vs_acc_between_joint, 
     SNR_vs_acc_stats, SNR_vs_acc_joint, SNR_vs_acc_joint_ci, SNR_vs_acc_joint_stats, 
     SNR_vs_session_stats, SNR_vs_session_joint, SNR_vs_session_joint_ci, SNR_vs_session_joint_stats,
     SNR_vs_accuracy_between_results, SNR_vs_accuracy_within_results, 
     SNR_vs_session_results,
     results, data_avg,
     file = file.path(r.path, 'multiverse_SNR.RData'))
