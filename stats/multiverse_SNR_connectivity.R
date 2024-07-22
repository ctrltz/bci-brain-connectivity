###
# Effects of SNR on connectivity
#
# Outputs:
# 1. fig6-snr-connectivity.png
# 2. fig6supp-connectivity-spectra.png
# 3. multiverse_SNR_connectivity.RData - intermediate results of the analysis
###

# Load the data
c(results, pipeline_desc) %<-% load_multiverse(input_filename)

# Create output folder for images
prefix = "multiverse_SNR_connectivity"
output.folder <- file.path(plot.path, prefix)
dir.create(output.folder)

# Set up the analysis
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

# Fit the models for all pipelines (split)
snr_vs_conn_stats <- lapply(param_space, fitMultiverseSplitLME, 
                            df = results, pipeline_desc = pipeline_desc)
snr_vs_conn_stats <- do.call(rbind, snr_vs_conn_stats)
snr_vs_conn_stats$Significant.MC <- snr_vs_conn_stats$p.value < p.mc
assert("PS ~ ROI SNR [split] did not converge",
       all(snr_vs_conn_stats$Converged))

# Fit a joint model
snr_vs_conn_joint_stats <- lapply(param_space, fitMultiverseJointLME, 
                                  df = results)
snr_vs_conn_joint_stats <- do.call(rbind, snr_vs_conn_joint_stats)
assert("PS ~ ROI SNR [joint] did not converge",
       all(snr_vs_conn_joint_stats$Converged))

# Combine the results
SNR_vs_conn_results <- lapply(param_space, getSummary,
                              split_stats = snr_vs_conn_stats,
                              joint_stats = snr_vs_conn_joint_stats,
                              commonFields = list(level = 'Within'))


### Multiverse Plot for SNR vs Connectivity --- Figure 6
snr_vs_conn_stats$fMeasure = factor(snr_vs_conn_stats$Measure, 
                                    levels = c('ImCoh', 'LagCoh', 'Coherence'))
snr_vs_conn_stats$fType = factor(snr_vs_conn_stats$Type, 
                                 levels = c('Within', 'Across'))
levels(snr_vs_conn_stats$fType) <- list("Within-hemisphere" = "Within",
                                        "Across-hemisphere" = "Across")
levels(snr_vs_conn_stats$Band) <- list("Broadband" = "BB", "Narrowband" = "NB")

p_snr_conn <- plotMultiverseSplit(snr_vs_conn_stats, 'Estimate', 'Significant.MC',
                                  facet_rule = 'fType + Band ~ fMeasure + Mask',
                                  val.name = bquote(beta))


### Sanity check: spectra of different connectivity measures
conn_df_no_mask <- conn_df[conn_df$Mask == "Anatomical",]
conn_df_mask <- conn_df[conn_df$Mask == "Task-based",]

# Anatomical ROIs -> main paper
p_imcoh <- plotConnectivitySpectra(conn_df_no_mask, "ImCoh", y.axis = TRUE)
p_lagcoh <- plotConnectivitySpectra(conn_df_no_mask, "LagCoh", x.axis = TRUE)
p_coh <- plotConnectivitySpectra(conn_df_no_mask, "Coherence")

# Task-based ROIs -> supplementary material
p_imcoh_mask <- plotConnectivitySpectra(conn_df_mask, "ImCoh", y.axis = TRUE)
p_lagcoh_mask <- plotConnectivitySpectra(conn_df_mask, "LagCoh", x.axis = TRUE)
p_coh_mask <- plotConnectivitySpectra(conn_df_mask, "Coherence")

### Combine the plots (figure 6)
legend <- get_legend(p_imcoh +
                       guides(linetype = guide_legend(title = "ROI Aggregation Method"),
                              color = guide_legend(title = 'Inverse Method')) +
                       theme(legend.direction = 'horizontal',
                             legend.position = 'bottom',
                             legend.key.width = unit(1, "cm")))
fig6A <- plot_grid(p_imcoh + theme(legend.position = 'none'),
                   p_lagcoh + theme(legend.position = 'none'),
                   p_coh + theme(legend.position = 'none'),
                   nrow = 1, rel_widths = c(0.35, 0.3, 0.35),
                   labels = c('A', '', ''))
fig6 <- plot_grid(fig6A, legend, p_snr_conn, 
                  ncol = 1, rel_heights = c(0.4, 0.05, 0.6),
                  labels = c('', '', 'B'))
save_plot(file.path(output.folder, 'fig6-snr-connectivity.pdf'),
          plot = fig6, base_width = 9, base_height = 11)
save_plot(file.path(output.folder, 'fig6-snr-connectivity.png'),
          plot = fig6, bg = "white", base_width = 9, base_height = 11)


### Combine the plots (supplementary)

legend <- get_legend(p_imcoh_mask +
                       guides(linetype = guide_legend(title = "ROI Aggregation Method"),
                              color = guide_legend(title = 'Inverse Method')) +
                       theme(legend.direction = 'horizontal',
                             legend.position = 'bottom',
                             legend.key.width = unit(1, "cm")))
fig6supp <- plot_grid(p_imcoh_mask + theme(legend.position = 'none'),
                      p_lagcoh_mask + theme(legend.position = 'none'),
                      p_coh_mask + theme(legend.position = 'none'),
                      nrow = 1, rel_widths = c(0.35, 0.3, 0.35))
fig6supp <- plot_grid(fig6supp, legend, 
                      ncol = 1, rel_heights = c(0.4, 0.05),
                      labels = c('', ''))
save_plot(file.path(output.folder, 'fig6supp-connectivity-spectra.pdf'),
          plot = fig6supp, base_width = 9, base_height = 5)
save_plot(file.path(output.folder, 'fig6supp-connectivity-spectra.png'),
          plot = fig6supp, bg = "white", base_width = 9, base_height = 6)


### Plot SNR vs Connectivity for different pipelines

for (ps in param_space) {
  message(paste('plotMultiverseScatterWithin:', ps$cols.scale[[1]], '~', ps$cols.scale[[2]]))
  p_separate <- plotMultiverseScatterWithin(results, ps$cols.scale[[1]], ps$cols.scale[[2]])
  ggsave(file.path(output.folder, paste('snr_vs_', ps$cols.scale[[2]], '.png', sep = '')),
         plot = p_separate, width = 9, height = 6)
}

### Save all the results
save(snr_vs_conn_stats, snr_vs_conn_joint_stats, SNR_vs_conn_results,
     file = file.path(r.path, 'multiverse_SNR_connectivity.RData'))