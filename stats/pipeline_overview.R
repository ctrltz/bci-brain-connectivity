###
# Prepare the pipeline overview figure (Figure 2)
#
# Outputs:
# 1. fig1-study-overview.png
# 2. fig2-multiverse-overview.png
# 3. overview.RData - underlying data (mu power contrasts, demo multiverse 
#      data frames)
###

# Create output folder for plots
prefix = "overview"
output.folder <- file.path(plot.path, prefix)
dir.create(output.folder)


## Figure 1 - Study Overview

# Training structure
p_training <- ggdraw() +
  draw_image(file.path(asset.path, 'training-structure.png'),
             x = 0.05, width = 0.9)

# Trial structure
p_trial <- ggdraw() +
  draw_image(file.path(asset.path, 'trial-structure.png'),
             x = 0.05, width = 0.9)


### Mu Power Contrasts
# Load mu power contrasts
mu_contrast <- sanity$mu.power.diff
c(n_subjects, n_periods, n_sensors) %<-% dim(mu_contrast)

# Load one EEG dataset to import channel locations properly
EEG <- import_set(file.path(eeg.path, 'S1_Session_1.set'))
assert("Channel locations are complete", nrow(EEG$chan_info) == 60)

# Calculate t-values for each period and EEG channel
t.vals <- data.frame()
for (p in 1:n_periods) {
  for (ch in 1:n_sensors) {
    t.vals <- rbind(t.vals, data.frame(Period = p, 
                                       x = EEG$chan_info$x[ch], 
                                       y = EEG$chan_info$y[ch], 
                                       amplitude = t.test(mu_contrast[,p,ch])$statistic[['t']],
                                       label = EEG$chan_info$electrode[ch]))
  }
}
t.vals$Period = as.factor(t.vals$Period)
levels(t.vals$Period) <- c("rest", "target", "feedback")

# Mu power contrast between imaginary movements of left and right hands
p_contrast <- ggplot(t.vals,
                     aes(x = x,
                         y = y,
                         fill = amplitude,
                         z = amplitude,
                         label = label)) +
  geom_topo(grid_res = 200,
            chan_size = rel(0.15), 
            head_size = rel(0.5),
            color = 'black',
            linetype = 'solid',
            linewidth = rel(0.1)) + 
  facet_wrap(. ~ Period, nrow = 1) +
  scale_fill_distiller(palette = "RdBu") + 
  theme_void() + 
  theme(legend.title = element_text(angle = -90), 
        legend.title.align = 0.5,
        strip.text.x = element_text(size = 10, margin = margin(b = 5))) +
  guides(fill = guide_colourbar(title.position = 'right',
                                title = 't-statistic')) +
  coord_equal()


# Laplacian
p_laplacian <- ggdraw() +
  draw_image(file.path(asset.path, 'laplacian.png'),
             x = 0.075, width = 0.85)

# Estimation of SNR with FOOOF
demo_snr <- sanity$demo.snr[,,1]
demo_df <- with(demo_snr, data.frame(
  PSD = t(spec.orig), noise = t(spec.noise),
  freqs = t(fit.freqs)))
mu_bins <- which(demo_snr$mu.bins > 0)
area_df <- with(demo_snr, data.frame(
  PSD = spec.orig[mu_bins],
  noise = spec.noise[mu_bins],
  freqs = t(mu.freqs)
))
p_demo_snr <- ggplot(data = demo_df, aes(x = freqs)) +
  geom_line(aes(y = PSD, color = 'PSD'), linewidth = 1) + 
  geom_line(aes(y = noise, color = 'noise'), linewidth = 1) +
  geom_ribbon(data = area_df, aes(x = freqs, ymin = noise, ymax = PSD, fill = 'Periodic'), alpha = 0.5) +
  geom_area(data = area_df, aes(x = freqs, y = noise, fill = 'Aperiodic'), alpha = 0.5) +
  geom_vline(xintercept = min(demo_snr$mu.freqs), color = 'gray', linetype = 'dashed') +
  geom_vline(xintercept = max(demo_snr$mu.freqs), color = 'gray', linetype = 'dashed') +
  scale_color_manual(values = c('black', 'gray'), limits = c('PSD', 'noise'),
                     labels = c('PSD', '1/f fit'), guide = 'none') +
  scale_fill_manual(values = c('orange', 'lightblue'),
                    limits = c('Periodic', 'Aperiodic')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab('Frequency (Hz)') + 
  ylab('PSD') +
  xlim(3, 33) +
  expand_limits(y = 0) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.left = element_blank())
p_snr_label <- ggdraw() +
  draw_label("Signal-to-Noise Ratio (SNR)", fontface = 'bold', y = 0.75) +
  draw_label("averaged over C3- and C4-Laplace", size = 12, y = 0.25)

# Combine the Laplace SNR plot
p_laplace_snr <- plot_grid(p_laplacian, p_demo_snr, 
                           nrow = 1, labels = c('D', 'E'))


# Within- and across-hemisphere phase synchronization
p_within_edges <- ggdraw() +
  draw_label("Within-hemisphere", size = 10, y = 0.925) +
  draw_image(file.path(asset.path, 'within-hemisphere.png'), 
             y = 0.05, height = 0.8)
p_across_edges <- ggdraw() +
  draw_label("Across-hemisphere", size = 10, y = 0.925) +
  draw_image(file.path(asset.path, 'across-hemisphere.png'), 
             y = 0.05, height = 0.8)
p_conn_label <- ggdraw() +
  draw_label("Phase Synchronization (PS)", fontface = 'bold', y = 0.75) +
  draw_label("between sensorimotor areas", size = 12, y = 0.25)

p_within_across <- plot_grid(p_within_edges, NULL, p_across_edges,
                             nrow = 1, rel_widths = c(1, 0.05, 1),
                             labels = c('F', '', ''))


p_bci_label <- ggdraw() +
  draw_label("Longitudinal BCI Training", fontface = 'bold', y = 0.75) +
  draw_label("cursor control based on motor imagery", size = 12, y = 0.25)


fig1_ABC <- plot_grid(p_bci_label, NULL, p_training, p_trial, p_contrast, 
                      ncol = 1, labels = c('', '', 'A', 'B', 'C'),
                      rel_heights = c(0.25, 0.05, 1, 1, 1))

fig1_DEF <- plot_grid(p_snr_label, NULL, p_laplace_snr, 
                      NULL, p_conn_label, p_within_across, 
                      ncol = 1,
                      rel_heights = c(0.25, 0.05, 1.415, 0.05, 0.25, 1.285))

fig1 <- plot_grid(fig1_ABC, NULL, fig1_DEF,
                  nrow = 1, rel_widths = c(1, 0.05, 1))

save_plot(file.path(output.folder, 'fig1-study-overview.pdf'), 
          fig1, base_width = 9, base_height = 6)
save_plot(file.path(output.folder, 'fig1-study-overview.png'), 
          fig1, bg = "white", base_width = 9, base_height = 6)


## Figure 2 - Overview of the Multiverse

# Multiverse of methods for source space analysis
p_source_space <- ggdraw() +
  draw_image(file.path(asset.path, 'source-space-analysis.png'),
             x = 0.05, width = 0.9)

# Anatomical and task-based ROI definitions
p_roi_defs_label_main <- ggdraw() +
  draw_label("ROI Definition", fontface = 'bold') 

p_roi_defs_label_anatomical <- ggdraw() +
  draw_label("Anatomical", size = 12)

p_roi_defs_label_task_based <- ggdraw() +
  draw_label("Task-based", size = 12)

p_roi_defs_anatomical <- ggdraw() +
  draw_image(file.path(asset.path, 'roi-definition-anatomical.png'))

p_roi_defs_task_based <- ggdraw() +
  draw_image(file.path(asset.path, 'roi-definition-task-based.png'))

# Load colors from MATLAB
colors <- readMat(color_filename)$colors
cm17 <- colors[[1]]
color_idx <- colors[[2]]
roi_colors <- rgb(cm17[color_idx, 1], cm17[color_idx, 2], 
                  cm17[color_idx, 3], maxColorValue = 1)
cm17 <- rgb(cm17[,1], cm17[,2], cm17[,3], maxColorValue = 1)

# Some dummy data for creating a suitable legend
df_defs <- data.frame(colors = roi_colors,
  labels = c(
  'Left precentral gyrus', 'Left postcentral gyrus',
  'Right precentral gyrus', 'Right postcentral gyrus'),
  x = c(1, 2, 3, 4))
p_roi_defs_misc <- ggplot(df_defs) +
  geom_bar(aes(x=x, fill=labels)) +
  scale_fill_manual(values=df_defs$colors, limits=df_defs$labels)
legend_defs <- get_legend(
  p_roi_defs_misc +
    guides(fill = guide_legend(ncol=2)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10))
)

p_roi_defs_main <- plot_grid(p_roi_defs_label_anatomical, 
                             p_roi_defs_label_task_based,
                             p_roi_defs_anatomical,
                             p_roi_defs_task_based,
                             nrow = 2, ncol = 2,
                             rel_heights = c(0.25, 1))
p_roi_defs_full <- plot_grid(p_roi_defs_label_main, 
                             p_roi_defs_main, NULL,  
                             legend_defs, 
                             ncol = 1, rel_heights = c(0.1, 1, 0.05, 0.3))

# Exemplary ROI aggregation weights (1SVD and AVG-F)
p_roi_weights_label_main <- ggdraw() +
  draw_label("ROI Aggregation Weights", fontface = 'bold') 

p_roi_weights_label_avg_flip <- ggdraw() +
  draw_label("AVG-F", size = 12)

p_roi_weights_label_1svd <- ggdraw() +
  draw_label("1SVD", size = 12)

p_roi_weights_avg_flip <- ggdraw() +
  draw_image(file.path(asset.path, 'roi-weights-avg-flip.png'),
             x = 0.1, width = 0.8)

p_roi_weights_1svd <- ggdraw() +
  draw_image(file.path(asset.path, 'roi-weights-1svd.png'),
             x = 0.1, width = 0.8)

# Some dummy data for creating a suitable colorbar
df_weights <- data.frame(x = c(1, 2, 3),
                         value = c(-1, 0, 1))
p_roi_weights_misc <- ggplot(df_weights) +
  geom_point(aes(x = x, y = value, fill = value), shape = 16) +
  scale_fill_gradientn(colors = cm17, limits = c(-1, 1))
legend_weights <- get_legend(
  p_roi_weights_misc +
    guides(fill = guide_colorbar(direction = "horizontal",
                                 title = "",
                                 title.hjust = 0.5,
                                 title.position = "top")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10))
)

p_roi_weights_main <- plot_grid(p_roi_weights_label_avg_flip, 
                                p_roi_weights_label_1svd,
                                p_roi_weights_avg_flip,
                                p_roi_weights_1svd,
                                nrow = 2, ncol = 2,
                                rel_heights = c(0.25, 1))
p_roi_weights_full <- plot_grid(p_roi_weights_label_main,
                                p_roi_weights_main, NULL,
                                legend_weights, 
                                ncol = 1, rel_heights = c(0.1, 1, 0.05, 0.3))


fig2_BC <- plot_grid(p_roi_defs_full, p_roi_weights_full, 
                     nrow = 1, labels = c('B', 'C'))



# Visualization example for the results of the split multiverse analysis
df_split_demo <- data.frame(
  Inverse = c(rep('eLORETA', 6), rep('LCMV', 6)),
  Band = rep('Broadband', 12),
  Mask = rep(c(rep('Anatomical', 3), rep('Task-based', 3)), 2),
  ROI_Method = rep(c('1SVD', '3SVD', 'AVG-F'), 4),
  Estimate = c(-0.13, -0.25, -0.11, -0.09, -0.19, -0.07,
               0.12, 0.2, 0.14, 0.1, 0.24, 0.08),
  p.value = c(1, 0.01, 1, 1, 0.01, 1,
              1, 0.01, 1, 1, 0.01, 1),
  Significant = c(F, T, F, F, T, F,
                  F, T, F, F, T, F)
)
p_multiverse_split_demo <- plotMultiverseSplit(df_split_demo, 'Estimate', 'Significant',
                                               lim = 0.3, facet_rule = 'Band ~ Mask')


# Visualization example for the results of the joint multiverse analysis
df_joint_demo <- data.frame(
  Category = rep(c('Average', rep('Within-/Across-hemisphere', times = 3)), times = 2),
  DummyX = rep('X', times = 8),
  DummyY = rep('Y', times = 8),
  Question = c(
    rep('Research\nQuestion 1', times = 4),
    rep('Research\nQuestion 2', times = 4)
  ),
  Measure = rep(c('SNR', 'ImCoh', 'LagCoh', 'Coherence'), times = 2),
  Consistency = c(0.9, 0.5, 0.4, 0.1, 0.9, 0.8, 0.8, 0.2),
  Significant.MC = c(T, T, T, T, F, F, F, T)
)
df_joint_demo$Measure <- factor(df_joint_demo$Measure,
                                levels = c('SNR', 'ImCoh', 'LagCoh', 'Coherence'))

p_multiverse_joint_demo <- plotMultiverseJoint(
  df_joint_demo, 'Consistency', 'Significant.MC', lim = 1, 
  facet_rule = 'Question ~ Category + Measure')


p_multiverse_label_main <- ggdraw() +
  draw_label("Multiverse Analysis - Examples of Visualization", fontface = 'bold') 

p_multiverse_label_split <- ggdraw() +
  draw_label("Split Analysis", size = 12)

p_multiverse_label_joint <- ggdraw() +
  draw_label("Joint Analysis", size = 12)


p_multiverse_split <- plot_grid(NULL, p_multiverse_split_demo, NULL,
                                nrow = 1, rel_widths = c(0.1, 1, 0.1))
p_multiverse_split <- plot_grid(p_multiverse_label_split, p_multiverse_split,
                                ncol = 1, rel_heights = c(0.1, 1))

p_multiverse_joint <- plot_grid(NULL, p_multiverse_joint_demo, NULL,
                                nrow = 1, rel_widths = c(0.05, 1, 0.05))
p_multiverse_joint <- plot_grid(p_multiverse_label_joint, p_multiverse_joint,
                                ncol = 1, rel_heights = c(0.1, 1))

fig2_DE <- plot_grid(p_multiverse_split, p_multiverse_joint, nrow = 1, 
                     labels = c('D', 'E'))

fig2 <- plot_grid(p_source_space, NULL, fig2_BC, NULL, p_multiverse_label_main, fig2_DE,
                  ncol = 1, rel_heights = c(1, 0, 1.2, 0.05, 0.2, 1.1), 
                  labels = c('A', '', ''))
save_plot(file.path(output.folder, 'fig2-pipeline-overview.pdf'), 
          fig2, base_width = 9, base_height = 10)
save_plot(file.path(output.folder, 'fig2-pipeline-overview.png'), 
          fig2, bg = "white", base_width = 9, base_height = 10)


# Save the results
save(t.vals, df_split_demo, df_joint_demo,
     file = file.path(r.path, 'overview.RData'))