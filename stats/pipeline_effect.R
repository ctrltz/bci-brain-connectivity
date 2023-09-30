###
# Analyze effects of different processing methods on the estimated values of 
# SNR and phase synchronization
#
# Outputs:
# 1. fig8-pipeline-effects-highlights.png - density plots of estimated values 
#    for different pipelines
# 2. fig8supp-snr-lcmv-eloreta.png - comparison of SNR with LCMV vs eLoreta
# 3. pipeline_effects.tex - t-values for effects of processing methods on the
#    estimated values (Table 2)
# 4. pipeline_effects.RData - intermediate results of the analysis
###

# Prefix for TeX output
prefix = "pipelineEffects"

# Load the data
c(results, results_SNR, pipeline_desc) %<-% load_multiverse(input_filename)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)


### Set up the plot theme
plot_theme <- theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1, fill = NA))


### Plot raw values for different pipelines

# SNR
plotPipelineHistograms(
  results_SNR, list(x ='LogSNR', name = '_snr_values.png'), 
  nrow = 2, plot.path = plot.path, prefix = prefix,
  width = 7, height = 4)

# Connectivity
args_list = list(
  list(x = 'ImCoh_Within', name = '_imcoh_within_values.png'),
  list(x = 'LagCoh_Within', name = '_lagcoh_within_values.png'),
  list(x = 'Coh_Within', name = '_coh_within_values.png'),
  list(x = 'ImCoh_Across', name = '_imcoh_across_values.png'),
  list(x = 'LagCoh_Across', name = '_lagcoh_across_values.png'),
  list(x = 'Coh_Across', name = '_coh_across_values.png')
)

lapply(args_list, plotPipelineHistograms,
       data = results, plot.path = plot.path, prefix = prefix,
       nrow = 4, width = 9, height = 6)


### Effect of the pipeline on the estimated values


estimatePipelineEffects <- function(df, params) {
  # Estimate effects of processing methods on the estimated PS values
  #
  # Parameters:
  #   df - data frame
  #   params - list (label, value, SNRCorr)
  #     label - row name for the future data frame with results
  #     value - value of interest
  #     SNRCorr - whether to include SNR as a covariate
  #
  # Returns:
  #   list(label, t-values, p-values)
  
  formula <- connectivity_vs_processing_formula(params$value, 
                                                SNRCorr = params$SNRCorr)
  message(paste('estimatePipelineEffects:', formula))
  
  fm_pipeline <- scaleAndFitLME(df,
    formula = formula,
    columns_to_scale = c(params$value)
  )
  assert(paste(params$value, "~ Processing did not converge"),
         has_converged(fm_pipeline))

  t.vals <- t(lmerTest:::get_coefmat(fm_pipeline)[,"t value"])
  p.vals <- t(lmerTest:::get_coefmat(fm_pipeline)[,"Pr(>|t|)"])
  
  list(params$label, t.vals, p.vals)
}


# SNR
SNR_multiverse <- scaleAndFitLME(
  results_SNR, 
  formula = snr_vs_processing_formula, 
  columns_to_scale = c('LogSNR')
)
assert("SNR ~ Processing [split] did not converge",
       has_converged(SNR_multiverse))
summary(SNR_multiverse)

t_SNR <- t(lmerTest:::get_coefmat(SNR_multiverse)[,"t value"])
p_SNR <- t(lmerTest:::get_coefmat(SNR_multiverse)[,"Pr(>|t|)"])


# Phase Synchronization
param_space = list(
  list(label = "WH ImCoh", value = "ImCoh_Within", SNRCorr = F),
  list(label = "WH LagCoh", value = "LagCoh_Within", SNRCorr = F),
  list(label = "WH Coherence", value = "Coh_Within", SNRCorr = F),
  list(label = "AH ImCoh", value = "ImCoh_Across", SNRCorr = F),
  list(label = "AH LagCoh", value = "LagCoh_Across", SNRCorr = F),
  list(label = "AH Coherence", value = "Coh_Across", SNRCorr = F),
  list(label = "WH ImCoh $|$ SNR", value = "ImCoh_Within", SNRCorr = T),
  list(label = "WH LagCoh $|$ SNR", value = "LagCoh_Within", SNRCorr = T),
  list(label = "WH Coherence $|$ SNR", value = "Coh_Within", SNRCorr = T),
  list(label = "AH ImCoh $|$ SNR", value = "ImCoh_Across", SNRCorr = T),
  list(label = "AH LagCoh $|$ SNR", value = "LagCoh_Across", SNRCorr = T),
  list(label = "AH Coherence $|$ SNR", value = "Coh_Across", SNRCorr = T)
)

pipelineEffects <- lapply(param_space, estimatePipelineEffects, df = results)
pipelineEffects <- c(list(list("SNR", t_SNR, p_SNR)), pipelineEffects)


### Summary
tvals <- as.data.frame(rbindlist(lapply(pipelineEffects, 
                                        function(x) as.data.frame(x[[2]])), 
                                 use.names = T, fill = T))
row.names(tvals) <- lapply(pipelineEffects, function(x) x[[1]])
pvals <- as.data.frame(rbindlist(lapply(pipelineEffects, 
                                        function(x) as.data.frame(x[[3]])), 
                                 use.names = T, fill = T))
row.names(pvals) <- lapply(pipelineEffects, function(x) x[[1]])


# Significant p-values
sig <- pvals %>%
  mutate_all(~ . < p.threshold)

# Significant p-values after correction for multiple comparisons (to be indicated with stars)
sig.mc <- pvals %>%
  mutate_all(~ if_else(. < p.mc, '*', '', missing = ''))

# Reformat t-values in the table: use bold for significant ones, add stars for multiple comparisons
reformat <- function(i) {
  col <- cols[[i]]
  
  if_else(
    is.na(tvals_latex[,col]),
    '---',
    cell_spec(paste(tvals_latex[,col], sig.mc[,col], sep = ''), 
              format = 'latex', bold = sig[,col])
  )
}

### Export to LaTeX
tvals_latex <- tvals
cols <- names(tvals_latex)
tvals_latex <- tvals_latex %>%
  mutate_if(is.numeric, round, 2)
tvals_latex <- as.data.frame(lapply(seq_along(tvals_latex), reformat))
names(tvals_latex) <- names(tvals)
rownames(tvals_latex) <- rownames(tvals)

tvals_latex %>%
  dplyr::select(-c("(Intercept)", "LogSNR")) %>%
  rename("3SVD $|$ 1SVD" = "ROI_Method3SVD", "AVG-F $|$ 1SVD" = "ROI_MethodAVG-F",
         "LCMV $|$ eLORETA" = "InverseLCMV", "With $|$ Without" = "MaskMask",
         "NB $|$ BB" = "BandNB") %>%
  knitr::kable(format = "latex", row.names = T, booktabs = T, escape = F, align = 'c',
               linesep = c("\\addlinespace", "", "", "\\addlinespace", "", "")) %>%
  add_header_above(c("Value", "Inverse Method", "ROI Method", 
                     "ROI Method", "Source Mask", "Band"), line = F,
                   align = c('l', 'c', 'c', 'c', 'c', 'c')) %>%
  tex.save(output_filename, "Summary", ., prefix = prefix)


### Plot a selection of the observed effects
set_scales <- list(
  scale_color_manual(values = c("eLORETA" = "#1984c5", "LCMV" = "#c23728")),
  scale_linetype_manual(values = c("1SVD" = "solid", "3SVD" = "dashed", "AVG-F" = "dotted")),
  scale_linewidth("Frequency Band", range = c(0.5, 1), breaks = c(1, 2), labels = c("Narrow", "Broad"))
)

# Use only results for anatomical ROIs to make the plots clearer
# Mask did not have a big effect anyways 
results_SNR_plot = results_SNR[results_SNR$Mask == "No Mask",]
results_plot = results[results$Mask == "No Mask",]

# Inverse -> SNR 
p1 <- ggplot(data = results_SNR_plot) +
  geom_line(mapping = aes(x = LogSNR, group = Pipeline, color = Inverse, linetype = ROI_Method),
            stat = "density", linewidth = 1) + 
  xlab('SNR (dB)') + ylab('Density') +
  facet_wrap(. ~ Inverse, nrow = 2) +
  ylim(0, NA) +
  set_scales +
  plot_theme

# Frequency Band -> ImCoh Within Hemispheres
band.labs <- c("Broadband", "Narrowband")
names(band.labs) <- c("BB", "NB")

p2 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = ImCoh_Within, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('ImCoh Within Hemispheres') + ylab('Density') +
  facet_wrap(. ~ Band, nrow = 2, 
             labeller = labeller(Band = band.labs)) +
  ylim(0, NA) +
  set_scales +
  plot_theme

# ROI_Method -> Coh Within Hemispheres
p3 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = Coh_Within, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('Coherence Within Hemispheres') + ylab('Density') +
  facet_wrap(. ~ ROI_Method, nrow = 3, scales = 'free_y') +
  ylim(0, NA) +
  set_scales +
  plot_theme

# Inverse -> ImCoh_Across
p4 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = ImCoh_Across, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('ImCoh Across Hemispheres') + ylab('Density') +
  facet_wrap(. ~ Inverse, nrow = 2) +
  ylim(0, NA) +
  set_scales +
  plot_theme

# Frequency Band -> ImCoh_Across
p5 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = ImCoh_Across, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('ImCoh Across Hemispheres') + ylab('Density') +
  facet_wrap(. ~ Band, nrow = 2, 
             labeller = labeller(Band = band.labs)) +
  ylim(0, NA) +
  set_scales +
  plot_theme

# ROI_Method -> Coh_Across
p6 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = Coh_Across, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('Coherence Across Hemispheres') + ylab('Density') +
  facet_wrap(. ~ ROI_Method, nrow = 3, scales = 'free_y') +
  ylim(0, NA) +
  set_scales +
  plot_theme


### Combine the plots --- Figure 8
prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  p4 + theme(legend.position="none"),
  p5 + theme(legend.position="none"),
  p6 + theme(legend.position="none"),
  align = 'vh',
  axis = 'l',
  labels = 'AUTO',
  hjust = -1,
  nrow = 2
)

# extract a legend that is laid out horizontally
legend <- get_legend(
  p2 + 
    guides(color = guide_legend(nrow=2, byrow=TRUE),
           linetype = guide_legend(nrow=2, byrow=TRUE),
           linewidth = guide_legend(nrow=2, byrow=TRUE)) +
    theme(legend.position = "bottom") +
    labs(color = "Inverse Method", linetype = "ROI Aggregation Method",
         linewidth = "Frequency Band")
)

# Combine all panels
fig8 <- plot_grid(prow, legend,
                  ncol = 1, axis = "h", labels = c('', '', 'D'),
                  rel_heights = c(1, .15))
save_plot(file.path(plot.path, 'fig8-pipeline-effects-highlights.pdf'), 
          fig8, base_width = 9, base_height = 7)
save_plot(file.path(plot.path, 'fig8-pipeline-effects-highlights.png'), 
          fig8, bg = "white", base_width = 9, base_height = 7)


### Scatter plot - compare SNR with eLORETA and LCMV (Figure 8 - supplementary)
eLoreta_pipelines <- grep('eLoreta', results_SNR$Pipeline, fixed = T)
LCMV_pipelines <- grep('LCMV', results_SNR$Pipeline, fixed = T)
cmp_df <- data.frame(
  Pipeline_eLoreta = results_SNR$Pipeline[eLoreta_pipelines],
  SNR_eLoreta = results_SNR$LogSNR[eLoreta_pipelines],
  Pipeline_LCMV = results_SNR$Pipeline[LCMV_pipelines],
  SNR_LCMV = results_SNR$LogSNR[LCMV_pipelines],
  Mask = results_SNR$Mask[eLoreta_pipelines],
  ROI_Method = results_SNR$ROI_Method[eLoreta_pipelines]
)
fig8supp <- ggplot(cmp_df) +
  geom_abline(intercept = 0, slope = 1, color = 'blue', linetype = 'dashed') +
  geom_point(aes(x = SNR_eLoreta, y = SNR_LCMV), color = 'darkgray', size = 0.5) +
  coord_fixed(ratio = 1) +
  xlab('SNR - eLoreta (dB)') + ylab('SNR - LCMV (dB)') +
  facet_grid(Mask ~ ROI_Method, switch = 'y') + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        strip.placement = 'outside')
save_plot(file.path(plot.path, 'fig8supp-snr-lcmv-eloreta.pdf'), 
          fig8supp, base_width = 9, base_height = 6)
save_plot(file.path(plot.path, 'fig8supp-snr-lcmv-eloreta.png'), 
          fig8supp, bg = "white", base_width = 9, base_height = 6)


### Save all the results
save(SNR_multiverse, pipelineEffects,
     tvals, pvals, cmp_df, file = file.path(r.path, 'pipeline_effects.RData'))