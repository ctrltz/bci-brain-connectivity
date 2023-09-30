###
# Prepare the pipeline overview figure (Figure 2)
#
# Outputs:
# 1. fig2-pipeline-overview.png
###

# Multiverse of methods for source space analysis
p_source_space <- ggdraw() +
  draw_image(file.path(asset.path, 'source-space-analysis.png'),
             x = 0.05, width = 0.9)

# Anatomical and task-based ROI definitions
p_roi_defs <- ggdraw() +
  draw_image(file.path(asset.path, 'fig2b-roi-definition.png'),
             x = 0.1, width = 0.8)

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

p_roi_defs_full <- plot_grid(p_roi_defs, legend_defs, 
                             ncol = 1, rel_heights = c(1, 0.3))

# Exemplary ROI aggregation weights (1SVD and AVG-flip)
p_roi_weights <- ggdraw() +
  draw_image(file.path(asset.path, 'fig2c-roi-weights.png'),
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
                                 title = "arb. units",
                                 title.hjust = 0.5,
                                 title.position = "top")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10))
)

p_roi_weights_full <- plot_grid(p_roi_weights, legend_weights, 
                                ncol = 1, rel_heights = c(1, 0.3))


fig2_BC <- plot_grid(p_roi_defs_full, p_roi_weights_full, 
                     nrow = 1, labels = c('B', 'C'))

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
  ylab(expression(paste('PSD (', mu, 'V'^2, '/Hz)'))) +
  xlim(3, 33) +
  expand_limits(y = 0) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.left = element_blank())
p_snr_label <- ggdraw() +
  draw_label("SNR", fontface = 'bold', y = 0.65)

fig2D <- plot_grid(NULL, p_snr_label, p_demo_snr,
                   ncol = 1, rel_heights = c(0.05, 0.15, 1), labels = c('D', ''))

# Within- and across-hemisphere phase synchronization
p_conn_edges <- ggdraw() +
  draw_image(file.path(asset.path, 'fig3e-within-across-edges.png'),
             x = 0.05, y = 0.05, width = 0.9)
p_conn_label <- ggdraw() +
  draw_label("Phase Synchronization (PS)", fontface = 'bold', y = 0.75) +
  draw_label("Within-hemisphere", size = 12, x = 0.25, y = 0.25) +
  draw_label("Across-hemisphere", size = 12, x = 0.75, y = 0.25)

fig2E <- plot_grid(NULL, p_conn_label, p_conn_edges,
                   ncol = 1, rel_heights = c(0.05, 0.25, 1), labels = c('E', ''))

# Visualization example for the results of the multiverse analysis
df_demo <- data.frame(
  Inverse = c(rep('eLORETA', 6), rep('LCMV', 6)),
  Band = rep('Broadband', 12),
  Mask = rep(c(rep('No Mask', 3), rep('Mask', 3)), 2),
  ROI_Method = rep(c('1SVD', '3SVD', 'AVG-F'), 4),
  Estimate = c(-0.13, -0.25, -0.11, -0.09, -0.19, -0.07,
               0.12, 0.2, 0.14, 0.1, 0.24, 0.08),
  p.value = c(1, 0.01, 1, 1, 0.01, 1,
              1, 0.01, 1, 1, 0.01, 1),
  Significant = c(F, T, F, F, T, F,
                  F, T, F, F, T, F)
)
p_multiverse_demo <- plotMultiverseSplit(df_demo, 'Estimate', 'Significant',
                                         lim = 0.3, facet_rule = 'Band ~ Mask')

p_multiverse_label <- ggdraw() +
  draw_label("Multiverse Analysis", fontface = 'bold', y = 0.75) +
  draw_label("Split Analysis - Visualization Example", size = 12, y = 0.25)

fig2F <- plot_grid(NULL, p_multiverse_label, p_multiverse_demo,
                   ncol = 1, rel_heights = c(0.05, 0.25, 1), labels = c('F', ''))

fig2_DEF <- plot_grid(fig2D, fig2E, fig2F, 
                      nrow = 1, rel_widths = c(1, 1.6, 1.4),
                      labels = c('', '', ''))

fig2 <- plot_grid(p_source_space, NULL, fig2_BC, NULL, fig2_DEF,
                  ncol = 1, rel_heights = c(1, 0, 1.2, 0.05, 1.2), 
                  labels = c('A', '', ''))
save_plot(file.path(plot.path, 'fig2-pipeline-overview.pdf'), 
          fig2, base_width = 9, base_height = 10)
save_plot(file.path(plot.path, 'fig2-pipeline-overview.png'), 
          fig2, bg = "white", base_width = 9, base_height = 10)