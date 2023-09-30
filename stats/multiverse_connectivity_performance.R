###
# Effects of connectivity on accuracy
#
# Outputs:
# 1. fig6-multiverse-connectivity-performance-within.png
# 2. fig6supp1-multiverse-connectivity-performance-between-subject.png
# 3. multiverse_connectivity_performance.RData - intermediate results of the 
#    analysis
###

### Between-Subject Effect of Connectivity Within Hemispheres on Accuracy

param_space <- list(
  list(x = 'ImCoh_Within', y = 'Accuracy', z = 'LogSNR', group = 'Subject',
       predictor = 'Average WH ImCoh', response = 'Average Accuracy',
       measure.name = 'ImCoh', measure.type = 'Within-Hemisphere'),
  list(x = 'LagCoh_Within', y = 'Accuracy', z = 'LogSNR', group = 'Subject',
       predictor = 'Average WH LagCoh', response = 'Average Accuracy',
       measure.name = 'LagCoh', measure.type = 'Within-Hemisphere'),
  list(x = 'Coh_Within', y = 'Accuracy', z = 'LogSNR', group = 'Subject',
       predictor = 'Average WH Coh', response = 'Average Accuracy',
       measure.name = 'Coherence', measure.type = 'Within-Hemisphere'),
  list(x = 'ImCoh_Across', y = 'Accuracy', z = 'LogSNR', group = 'Subject',
       predictor = 'Average AH ImCoh', response = 'Average Accuracy',
       measure.name = 'ImCoh', measure.type = 'Across-Hemisphere'),
  list(x = 'LagCoh_Across', y = 'Accuracy', z = 'LogSNR', group = 'Subject',
       predictor = 'Average AH LagCoh', response = 'Average Accuracy',
       measure.name = 'LagCoh', measure.type = 'Across-Hemisphere'),
  list(x = 'Coh_Across', y = 'Accuracy', z = 'LogSNR', group = 'Subject',
       predictor = 'Average AH Coh', response = 'Average Accuracy',
       measure.name = 'Coherence', measure.type = 'Across-Hemisphere')
)

conn_vs_acc_between_nosnr_stats <- lapply(param_space, fitMultiverseSplitBetween,
                                          df = results, pipeline_desc = pipeline_desc)
conn_vs_acc_between_nosnr_stats <- do.call(rbind, conn_vs_acc_between_nosnr_stats)
conn_vs_acc_between_nosnr_stats$Correction <- 'Not Corrected for SNR'
conn_vs_acc_between_nosnr_stats$Significant.MC <- conn_vs_acc_between_nosnr_stats$p.value < p.mc


conn_vs_acc_between_snr_stats <- lapply(param_space, fitMultiverseSplitBetweenPartial,
                                        df = results, pipeline_desc = pipeline_desc)
conn_vs_acc_between_snr_stats <- do.call(rbind, conn_vs_acc_between_snr_stats)
conn_vs_acc_between_snr_stats$Correction <- 'Corrected for SNR'
conn_vs_acc_between_snr_stats$Significant.MC <- conn_vs_acc_between_snr_stats$p.value < p.mc

conn_vs_acc_between_stats <- rbind(conn_vs_acc_between_nosnr_stats,
                                   conn_vs_acc_between_snr_stats)

# Plot the multiverse
lim_conn_acc_between = ceilingn(max(abs(conn_vs_acc_between_stats$Estimate)), 2)

conn_vs_acc_between_stats <- within(conn_vs_acc_between_stats, {
  fMeasure = factor(Measure, levels = c('ImCoh', 'LagCoh', 'Coherence'))
  Band = factor(Band, levels = c('BB', 'NB'), labels = c('Broadband', 'Narrowband'))
  fType = factor(Type)
  fCorrection = factor(Correction, levels = c('Not Corrected for SNR', 'Corrected for SNR'))
})

p_conn_acc_between_within <- plotMultiverseSplit(
  conn_vs_acc_between_stats[conn_vs_acc_between_stats$fType == 'Within-Hemisphere',],
  'Estimate', 'Significant.MC', facet_rule = 'fCorrection + Band ~ fMeasure + Mask',
  val.name = '(Partial) Correlation', lim = lim_conn_acc_between)
p_conn_acc_between_across <- plotMultiverseSplit(
  conn_vs_acc_between_stats[conn_vs_acc_between_stats$fType == 'Across-Hemisphere',],
  'Estimate', 'Significant.MC', facet_rule = 'fCorrection + Band ~ fMeasure + Mask',
  val.name = '(Partial) Correlation', lim = lim_conn_acc_between)

### Combine the figures
fig6supp1 <- plot_grid(
  p_conn_acc_between_within + theme(legend.position = 'none'),
  p_conn_acc_between_across + theme(legend.position = 'none'),
  get_legend(p_conn_acc_between_within),
  ncol = 1, rel_heights = c(1, 1, 0.2),
  labels = c('A', 'B', '')
)

save_plot(file.path(plot.path, prefix, 'fig6supp1-multiverse-connectivity-performance-between-subject.png'),
          fig6supp1, bg = "white", base_width = 8, base_height = 12)



### Within-Subject Effect of Connectivity Within Hemispheres on Accuracy

param_space <- list(
  list(formula_split = accuracy_vs_connectivity_formula("ImCoh_Within"),
       formula_joint = accuracy_vs_connectivity_formula("ImCoh_Within", joint = T),
       cols.scale = c('Accuracy', 'ImCoh_Within'), 
       predictor = 'ImCoh_Within', response = 'Accuracy', 
       measure.name = 'ImCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("LagCoh_Within"),
       formula_joint = accuracy_vs_connectivity_formula("LagCoh_Within", joint = T),
       cols.scale = c('Accuracy', 'LagCoh_Within'), 
       predictor = 'LagCoh_Within', response = 'Accuracy', 
       measure.name = 'LagCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("Coh_Within"),
       formula_joint = accuracy_vs_connectivity_formula("Coh_Within", joint = T),
       cols.scale = c('Accuracy', 'Coh_Within'), 
       predictor = 'Coh_Within', response = 'Accuracy', 
       measure.name = 'Coherence', measure.type = 'Not Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("ImCoh_Within", SNRCorr = T),
       formula_joint = accuracy_vs_connectivity_formula("ImCoh_Within", SNRCorr = T, joint = T),
       cols.scale = c('Accuracy', 'LogSNR', 'ImCoh_Within'), 
       predictor = 'ImCoh_Within', response = 'Accuracy $|$ SNR',
       measure.name = 'ImCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("LagCoh_Within", SNRCorr = T),
       formula_joint = accuracy_vs_connectivity_formula("LagCoh_Within", SNRCorr = T, joint = T),
       cols.scale = c('Accuracy', 'LogSNR', 'LagCoh_Within'), 
       predictor = 'LagCoh_Within', response = 'Accuracy $|$ SNR', 
       measure.name = 'LagCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("Coh_Within", SNRCorr = T),
       formula_joint = accuracy_vs_connectivity_formula("Coh_Within", SNRCorr = T, joint = T),
       cols.scale = c('Accuracy', 'LogSNR', 'Coh_Within'), 
       predictor = 'Coh_Within', response = 'Accuracy $|$ SNR',
       measure.name = 'Coherence', measure.type = 'Corrected for SNR')
)

conn_within_vs_acc_stats <- lapply(param_space, fitMultiverseSplitLME,
                                   df = results, pipeline_desc = pipeline_desc)
conn_within_vs_acc_stats <- do.call(rbind, conn_within_vs_acc_stats)
conn_within_vs_acc_stats$Significant.MC <- conn_within_vs_acc_stats$p.value < p.mc
assert("Accuracy ~ Within PS [split] did not converge",
       all(conn_within_vs_acc_stats$Converged))

conn_within_vs_acc_joint_stats <- lapply(param_space, fitMultiverseJointLME,
                                         df = results)
conn_within_vs_acc_joint_stats <- do.call(rbind, conn_within_vs_acc_joint_stats)
assert("Accuracy ~ Within PS [joint] did not converge",
       all(conn_within_vs_acc_joint_stats$Converged))

conn_within_vs_acc_results <- lapply(param_space, getSummary,
                                     split_stats = conn_within_vs_acc_stats,
                                     joint_stats = conn_within_vs_acc_joint_stats,
                                     commonFields = list(type2 = 'Within', 
                                                         level = 'Within'))
conn_within_vs_acc_results <- lapply(conn_within_vs_acc_results,
                                     \(x) {
                                       x$SNRCorr <- x$type
                                       x$type <- x$type2
                                       x$type2 <- NULL
                                       x
                                     })


### ImCoh Within Hemispheres vs Accuracy
p_imcoh_within_acc <- ggplot(data = results) +
  geom_point(mapping = aes(x = ImCoh_Within, y = Accuracy, group = Subject, color = Subject), size = 0.5) +
  geom_line(mapping = aes(x = ImCoh_Within, y = Accuracy, group = Subject), color = "black", alpha = 0.5, 
            stat="smooth", method = "lm", se = F, linewidth = 0.5) +
  geom_line(mapping = aes(x = ImCoh_Within, y = Accuracy), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 4) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")

ggsave(file.path(plot.path, prefix, 'multiverse_imcoh_within_vs_accuracy.png'),
       p_imcoh_within_acc, width = 9, height = 6)


### Coherence Within Hemispheres vs Accuracy
p_coh_within_acc <- ggplot(data = results) +
  geom_point(mapping = aes(x = Coh_Within, y = Accuracy, group = Subject, color = Subject), size = 0.5) +
  geom_line(mapping = aes(x = Coh_Within, y = Accuracy, group = Subject), color = "black", alpha = 0.5, 
            stat="smooth", method = "lm", se = F, linewidth = 0.5) +
  geom_line(mapping = aes(x = Coh_Within, y = Accuracy), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 4) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")

ggsave(file.path(plot.path, prefix, 'multiverse_coh_within_vs_accuracy.png'),
       p_coh_within_acc, width = 9, height = 6)


### Within-Subject Effect of Connectivity Across Hemispheres on Accuracy

param_space <- list(
  list(formula_split = accuracy_vs_connectivity_formula("ImCoh_Across"),
       formula_joint = accuracy_vs_connectivity_formula("ImCoh_Across", joint = T),
       cols.scale = c('Accuracy', 'ImCoh_Across'), 
       predictor = 'ImCoh_Across', response = 'Accuracy',
       measure.name = 'ImCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("LagCoh_Across"),
       formula_joint = accuracy_vs_connectivity_formula("LagCoh_Across", joint = T),
       cols.scale = c('Accuracy', 'LagCoh_Across'), 
       predictor = 'LagCoh_Across', response = 'Accuracy',
       measure.name = 'LagCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("Coh_Across"),
       formula_joint = accuracy_vs_connectivity_formula("Coh_Across", joint = T),
       cols.scale = c('Accuracy', 'Coh_Across'), 
       predictor = 'Coh_Across', response = 'Accuracy',
       measure.name = 'Coherence', measure.type = 'Not Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("ImCoh_Across", SNRCorr = T),
       formula_joint = accuracy_vs_connectivity_formula("ImCoh_Across", SNRCorr = T, joint = T),
       cols.scale = c('Accuracy', 'LogSNR', 'ImCoh_Across'), 
       predictor = 'ImCoh_Across', response = 'Accuracy $|$ SNR',
       measure.name = 'ImCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("LagCoh_Across", SNRCorr = T),
       formula_joint = accuracy_vs_connectivity_formula("LagCoh_Across", SNRCorr = T, joint = T),
       cols.scale = c('Accuracy', 'LogSNR', 'LagCoh_Across'), 
       predictor = 'LagCoh_Across', response = 'Accuracy $|$ SNR',
       measure.name = 'LagCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = accuracy_vs_connectivity_formula("Coh_Across", SNRCorr = T),
       formula_joint = accuracy_vs_connectivity_formula("Coh_Across", SNRCorr = T, joint = T),
       cols.scale = c('Accuracy', 'LogSNR', 'Coh_Across'),
       predictor = 'Coh_Across', response = 'Accuracy $|$ SNR',
       measure.name = 'Coherence', measure.type = 'Corrected for SNR')
)

conn_across_vs_acc_stats <- lapply(param_space, fitMultiverseSplitLME,
                                   df = results, pipeline_desc = pipeline_desc)
conn_across_vs_acc_stats <- do.call(rbind, conn_across_vs_acc_stats)
conn_across_vs_acc_stats$Significant.MC <- conn_across_vs_acc_stats$p.value < p.mc
assert("Accuracy ~ Across PS [split] did not converge",
       all(conn_across_vs_acc_stats$Converged))

conn_across_vs_acc_joint_stats <- lapply(param_space, fitMultiverseJointLME,
                                         df = results)
conn_across_vs_acc_joint_stats <- do.call(rbind, conn_across_vs_acc_joint_stats)
assert("Accuracy ~ Across PS [joint] did not converge",
       all(conn_across_vs_acc_joint_stats$Converged))

conn_across_vs_acc_results <- lapply(param_space, getSummary,
                                     split_stats = conn_across_vs_acc_stats,
                                     joint_stats = conn_across_vs_acc_joint_stats,
                                     commonFields = list(type2 = 'Across', 
                                                         level = 'Within'))
conn_across_vs_acc_results <- lapply(conn_across_vs_acc_results,
                                     \(x) {
                                       x$SNRCorr <- x$type
                                       x$type <- x$type2
                                       x$type2 <- NULL
                                       x
                                     })



### Coherence Across Hemispheres vs Accuracy
p_coh_across_acc <- ggplot(data = results) +
  geom_point(mapping = aes(x = Coh_Across, y = Accuracy, group = Subject, color = Subject), size = 0.5) +
  geom_line(mapping = aes(x = Coh_Across, y = Accuracy, group = Subject), color = "black", alpha = 0.5, 
            stat="smooth", method = "lm", se = F, linewidth = 0.5) +
  geom_line(mapping = aes(x = Coh_Across, y = Accuracy), color = "blue",
            stat="smooth", method = "lm", se = F) +
  facet_wrap(. ~ Pipeline, nrow = 4) + 
  scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")

ggsave(file.path(plot.path, prefix, 'multiverse_coh_across_vs_accuracy.png'),
       p_coh_across_acc, width = 9, height = 6)


### Plot the Multiverse
lim_conn_acc_within = ceilingn(max(abs(c(conn_within_vs_acc_stats$Estimate,
                                         conn_across_vs_acc_stats$Estimate))), 2)

# Within-hemisphere connectivity
conn_within_vs_acc_stats <- within(conn_within_vs_acc_stats, {
  fMeasure = factor(Measure, levels = c('ImCoh', 'LagCoh', 'Coherence'))
  Band = factor(Band, levels = c('BB', 'NB'), labels = c('Broadband', 'Narrowband'))
  fType = factor(Type, levels = c('Not Corrected for SNR', 'Corrected for SNR'))
})

p_conn_acc_within <- plotMultiverseSplit(conn_within_vs_acc_stats, 'Estimate', 
                                         'Significant.MC', lim = lim_conn_acc_within, 
                                         facet_rule = 'fType + Band ~ fMeasure + Mask')

# Across-hemisphere connectivity
conn_across_vs_acc_stats <- within(conn_across_vs_acc_stats, {
  fMeasure = factor(Measure, levels = c('ImCoh', 'LagCoh', 'Coherence'))
  Band = factor(Band, levels = c('BB', 'NB'), labels = c('Broadband', 'Narrowband'))
  fType = factor(Type, levels = c('Not Corrected for SNR', 'Corrected for SNR'))
})

p_conn_acc_across <- plotMultiverseSplit(conn_across_vs_acc_stats, 'Estimate', 
                                         'Significant.MC', lim = lim_conn_acc_within, 
                                         facet_rule = 'fType + Band ~ fMeasure + Mask')


### Combine the figures
fig6 <- plot_grid(
  p_conn_acc_within + theme(legend.position = 'none'),
  p_conn_acc_across + theme(legend.position = 'none'),
  get_legend(p_conn_acc_within),
  ncol = 1, rel_heights = c(1, 1, 0.2),
  labels = c('A', 'B', '')
)

save_plot(file.path(plot.path, prefix, 'fig6-multiverse-connectivity-performance-within.png'),
          fig6, bg = "white", base_width = 8, base_height = 12)


### Save all the results
save(conn_vs_acc_between_stats, conn_within_vs_acc_stats, conn_within_vs_acc_joint_stats,
     conn_across_vs_acc_stats, conn_across_vs_acc_joint_stats,
     conn_within_vs_acc_results, conn_across_vs_acc_results,
     file = file.path(r.path, 'multiverse_connectivity_performance.RData'))