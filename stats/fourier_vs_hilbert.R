### 
# Compare Fourier-based vs Hilbert-based PS for narrowband pipelines to show
# systematic differences
#
# Outputs:
# 1. figR1-fourier-vs-hilbert.png - scatter plots comparing the estimates for
#    all narrowband pipelines (only included in the response letter)
###

# Create output folder
output.folder <- file.path(plot.path, 'revision')
dir.create(output.folder)

# Load data for Fourier-only estimation vs Fourier and Hilbert for comparison
bb_only_filename <- file.path(data.path, 'BCI_MI_multiverse_rest_results_long.mat')
bb_nb_filename <- file.path(data.path, 'BCI_MI_multiverse_rest_results_long_BB_NB.mat')
  
c(results_bb_only, pipeline_desc) %<-% load_multiverse(bb_only_filename)
c(results_bb_nb, pipeline_desc) %<-% load_multiverse(bb_nb_filename)

# Sanity check: BB should show no difference
# results_bb_only_nb <- results_bb_only[grepl('bb', results_bb_only$Pipeline, fixed = T),]
# results_bb_nb_nb <- results_bb_nb[grepl('bb', results_bb_nb$Pipeline, fixed = T),]

results_bb_only_nb <- results_bb_only[grepl('nb', results_bb_only$Pipeline, fixed = T),]
results_bb_nb_nb <- results_bb_nb[grepl('nb', results_bb_nb$Pipeline, fixed = T),]
results_cmp <- merge(results_bb_only_nb, results_bb_nb_nb, by = 
  c('Pipeline', 'Subject', 'Session', 'Accuracy', 'Inverse', 'Band', 
    'ROI_Agg', 'Mask', 'ROI_Method', 'BandViz'), suffixes = c('_BB', '_NB'))


p_imcoh_within <- ggplot(results_cmp) +
  geom_point(aes(x = ImCoh_Within_BB, y = ImCoh_Within_NB), size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  facet_nested(ROI_Method ~ Inverse + Mask) +
  xlab('Fourier-based estimate') + ylab('Hilbert-based estimate') +
  ggtitle('Within-hemisphere ImCoh') +
  coord_fixed() +
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


p_imcoh_across <- ggplot(results_cmp) +
  geom_point(aes(x = ImCoh_Across_BB, y = ImCoh_Across_NB), size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  facet_nested(ROI_Method ~ Inverse + Mask) +
  xlab('Fourier-based estimate') + ylab('Hilbert-based estimate') +
  expand_limits(x = 0) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  ggtitle('Across-hemisphere ImCoh') +
  coord_fixed() +
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


p_lagcoh_within <- ggplot(results_cmp) +
  geom_point(aes(x = LagCoh_Within_BB, y = LagCoh_Within_NB), size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  facet_nested(ROI_Method ~ Inverse + Mask) +
  xlab('Fourier-based estimate') + ylab('Hilbert-based estimate') +
  expand_limits(x = 0) +
  ggtitle('Within-hemisphere LagCoh') +
  coord_fixed() +
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


p_lagcoh_across <- ggplot(results_cmp) +
  geom_point(aes(x = LagCoh_Across_BB, y = LagCoh_Across_NB), size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  facet_nested(ROI_Method ~ Inverse + Mask) +
  xlab('Fourier-based estimate') + ylab('Hilbert-based estimate') +
  expand_limits(x = 0) +
  ggtitle('Across-hemisphere LagCoh') +
  coord_fixed() +
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


p_coh_within <- ggplot(results_cmp) +
  geom_point(aes(x = Coh_Within_BB, y = Coh_Within_NB), size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  facet_nested(ROI_Method ~ Inverse + Mask) +
  xlab('Fourier-based estimate') + ylab('Hilbert-based estimate') +
  expand_limits(x = 0) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  ggtitle('Within-hemisphere Coherence') +
  coord_fixed() +
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


p_coh_across <- ggplot(results_cmp) +
  geom_point(aes(x = Coh_Across_BB, y = Coh_Across_NB), size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  facet_nested(ROI_Method ~ Inverse + Mask) +
  xlab('Fourier-based estimate') + ylab('Hilbert-based estimate') +
  expand_limits(x = 0) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  ggtitle('Across-hemisphere Coherence') +
  coord_fixed() +
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


p_fourier_vs_hilbert <- plot_grid(p_imcoh_within, p_imcoh_across,
                                  p_lagcoh_within, p_lagcoh_across,
                                  p_coh_within, p_coh_across,
                                  ncol = 2, labels = 'AUTO')
save_plot(file.path(output.folder, 'figR1-fourier-vs-hilbert.pdf'), 
          p_fourier_vs_hilbert, base_width = 9, base_height = 12)
save_plot(file.path(output.folder, 'figR1-fourier-vs-hilbert.png'), 
          p_fourier_vs_hilbert, bg = "white", base_width = 9, base_height = 12)