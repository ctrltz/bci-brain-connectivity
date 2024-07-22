###
# Analyze effects of different processing methods on the estimated values of 
# SNR and phase synchronization
#
# Outputs:
# 1. fig9-pipeline-effects-highlights.png - density plots of estimated values 
#    for different pipelines
# 2. fig9supp1-snr-band-inverse-interaction.png - interaction between filtering 
#      and inverse modeling
# 3. fig9supp2-coh-within-mask-roi-method-interaction.png - interaction between 
#      ROI definition and ROI aggregation methods
# 4. fig9supp3-snr-inverse-band-interaction.png - another way to plot the
#      interaction between filtering and inverse modeling
# 5. pipeline_effects.tex - t-values for effects of processing methods on the
#    estimated values (Tab. 2) and list of significant interactions (Tab. S6)
# 6. pipeline_effects.RData - intermediate results of the analysis
###

# Prefix for TeX output
prefix = "pipelineEffects"

# Create output folder for images
output.folder <- file.path(plot.path, "pipeline_effects")
dir.create(output.folder)

# Load the data
c(results, pipeline_desc) %<-% load_multiverse(input_filename)

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
  results, list(x ='LogSNR', name = '_snr_values.png'), 
  nrow = 4, plot.path = output.folder, prefix = "pipeline_effects",
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
       data = results, plot.path = output.folder, prefix = "pipeline_effects",
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
message(paste('estimatePipelineEffects:', snr_vs_processing_formula))
SNR_multiverse <- scaleAndFitLME(
  results, 
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


# Split into main effects and interactions
main_idx <- c(1:6, 16)
int_idx <- 7:15

tvals_main <- tvals[,main_idx]
tvals_int <- tvals[,int_idx]
pvals_main <- pvals[,main_idx]
pvals_int <- pvals[,int_idx]


# Report main effects

# Significant p-values
sig_main <- pvals_main %>%
  mutate_all(~ . < p.threshold)

# Significant p-values after correction for multiple comparisons (to be indicated with stars)
sig_main.mc <- pvals_main %>%
  mutate_all(~ if_else(. < p.mc, '*', '', missing = ''))

# Reformat t-values in the table: use bold for significant ones, add stars for multiple comparisons
reformat <- function(i, df, df_sig, df_sig_mc) {
  col <- cols[[i]]
  
  if_else(
    is.na(df[,col]),
    '---',
    cell_spec(paste(df[,col], df_sig_mc[,col], sep = ''), 
              format = 'latex', bold = df_sig[,col])
  )
}

### Export to LaTeX
tvals_latex <- tvals_main
cols <- names(tvals_latex)
tvals_latex <- tvals_latex %>%
  mutate_if(is.numeric, round, 2)
tvals_latex <- as.data.frame(lapply(seq_along(tvals_latex), reformat, 
                                    tvals_latex, sig_main, sig_main.mc))
names(tvals_latex) <- names(tvals_main)
rownames(tvals_latex) <- rownames(tvals_main)

tex.save(output_filename, "% Main effects\n")
tvals_latex %>%
  dplyr::select(-c("(Intercept)", "LogSNR")) %>%
  rename("3SVD $|$ 1SVD" = "ROI_Method3SVD", "AVG-F $|$ 1SVD" = "ROI_MethodAVG-F",
         "LCMV $|$ eLORETA" = "InverseLCMV", "Task $|$ Anat." = "MaskTask-based",
         "NB $|$ BB" = "BandNB") %>%
  knitr::kable(format = "latex", row.names = T, booktabs = T, escape = F, align = 'c',
               linesep = c("\\addlinespace", "", "", "\\addlinespace", "", "")) %>%
  add_header_above(c("Value", "Band", "Inverse Method", 
                     "ROI Def.", "ROI Method", "ROI Method"), line = F,
                   align = c('l', 'c', 'c', 'c', 'c', 'c')) %>%
  tex.save(output_filename, "Summary", ., prefix = prefix)


# Report the interaction effects

tvals_int_wide <- cbind("Response" = rownames(tvals_int), tvals_int)
tvals_int_wide <- tvals_int_wide %>% pivot_longer(
  cols = !Response, names_to="Interaction", values_to="t-value")

pvals_int_wide <- cbind("Response" = rownames(pvals_int), pvals_int)
pvals_int_wide <- pvals_int_wide %>% pivot_longer(
  cols = !Response, names_to="Interaction", values_to="p-value")

stats_int_wide <- tvals_int_wide
stats_int_wide$`p-value` <- pvals_int_wide$`p-value`
stats_int_wide$sig <- stats_int_wide$`p-value` < p.threshold
stats_int_wide$sig.mc <- if_else(stats_int_wide[,'p-value'] < p.mc, '*', '', missing = '')

stats_int_latex <- stats_int_wide %>%
  arrange(Interaction) %>%
  filter(`p-value` < p.threshold) %>%
  mutate_at(vars(`t-value`, `p-value`), round, 2)
stats_int_latex$Interaction <- with(stats_int_latex, case_when(
  Interaction == "BandNB:InverseLCMV" ~ "Band (NB) $\\cdot$ Inverse (LCMV)",
  Interaction == "InverseLCMV:MaskTask-based" ~ "Inverse (LCMV) $\\cdot$ ROI Definition (Task-based)",
  Interaction == "MaskTask-based:ROI_Method3SVD" ~ "ROI Definition (Task-based) $\\cdot$ ROI Method (3SVD)",
  Interaction == "MaskTask-based:ROI_MethodAVG-F" ~ "ROI Definition (Task-based) $\\cdot$ ROI Method (AVG-F)",
  .default = Interaction
))
stats_int_latex$`p-value` <- with(stats_int_latex, cell_spec(
  paste(if_else(`p-value` < 0.001, '$<$0.001', 
                format(`p-value`, digits = 3)), 
        sig.mc, sep = ''), 
  format = 'latex', bold = sig, escape = F))

tex.save(output_filename, "\n% Interactions\n")
stats_int_latex %>%
  dplyr::select(Interaction, Response, `t-value`, `p-value`) %>%
  knitr::kable(format = "latex", row.names = , booktabs = T, escape = F, align = 'c') %>%
  collapse_rows(column = 1) %>%
  gsub("\\cmidrule{2-4}\n", "", ., fixed = T) %>%
  tex.save(output_filename, "InteractionSummary", ., prefix = prefix)


  
### Plot a selection of the observed effects
set_scales <- list(
  scale_color_manual(values = c("eLORETA" = "#1984c5", "LCMV" = "#c23728")),
  scale_linetype_manual(values = c("1SVD" = "solid", "3SVD" = "dashed", "AVG-F" = "dotted")),
  scale_linewidth("Frequency Band", range = c(0.5, 1), breaks = c(1, 2), labels = c("Narrow", "Broad"))
)

# Use only results for anatomical ROIs to make the plots clearer
results_plot = results[results$Mask == "Anatomical",]

# Inverse -> SNR 
p1 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = LogSNR, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
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
  xlab('Within-hemisphere ImCoh') + ylab('Density') +
  facet_wrap(. ~ ROI_Method, nrow = 3, scales = 'free_y') +
  ylim(0, NA) +
  set_scales +
  plot_theme

# ROI_Method -> Coh Within Hemispheres
p3 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = Coh_Within, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('Within-hemisphere Coherence') + ylab('Density') +
  facet_wrap(. ~ ROI_Method, nrow = 3, scales = 'free_y') +
  ylim(0, NA) +
  set_scales +
  plot_theme

# Inverse -> Coh_Across
p4 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = Coh_Across, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('Across-hemisphere Coherence') + ylab('Density') +
  facet_wrap(. ~ Inverse, nrow = 2) +
  ylim(0, NA) +
  set_scales +
  plot_theme

# Frequency Band -> ImCoh_Across
p5 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = ImCoh_Across, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('Across-hemisphere ImCoh') + ylab('Density') +
  facet_wrap(. ~ ROI_Method, nrow = 3, scales = 'free_y') +
  ylim(0, NA) +
  set_scales +
  plot_theme

# ROI_Method -> Coh_Across
p6 <- ggplot(data = results_plot) +
  geom_line(mapping = aes(x = Coh_Across, group = Pipeline, 
                          color = Inverse, linewidth = BandViz, linetype = ROI_Method),
            stat = "density") + 
  xlab('Across-hemisphere Coherence') + ylab('Density') +
  facet_wrap(. ~ ROI_Method, nrow = 3, scales = 'free_y') +
  ylim(0, NA) +
  set_scales +
  plot_theme


### Combine the plots --- Figure 9
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
fig9 <- plot_grid(prow, legend,
                  ncol = 1, axis = "h", labels = c('', '', 'D'),
                  rel_heights = c(1, .15))
save_plot(file.path(output.folder, 'fig9-pipeline-effects-highlights.pdf'), 
          fig9, base_width = 9, base_height = 7)
save_plot(file.path(output.folder, 'fig9-pipeline-effects-highlights.png'), 
          fig9, bg = "white", base_width = 9, base_height = 7)


### Plot all significant interactions
stats_int_plot <- stats_int_wide %>%
  mutate(value_col = recode(Response, "SNR" = "LogSNR", "WH ImCoh" = "ImCoh_Within",
                            "WH LagCoh" = "LagCoh_Within", "WH Coherence" = "Coh_Within",
                            "AH ImCoh" = "ImCoh_Across", "AH LagCoh" = "LagCoh_Across",
                            "AH Coherence" = "Coh_Across", .default = "")) %>%
  mutate(group_col = recode(Interaction,
                            "BandNB:InverseLCMV" = "Inverse",
                            "InverseLCMV:MaskTask-based" = "Inverse",
                            "MaskTask-based:ROI_Method3SVD" = "Mask",
                            "MaskTask-based:ROI_MethodAVG-F" = "Mask")) %>%
  mutate(other_col = recode(Interaction,
                            "BandNB:InverseLCMV" = "Band",
                            "InverseLCMV:MaskTask-based" = "Mask",
                            "MaskTask-based:ROI_Method3SVD" = "ROI_Method",
                            "MaskTask-based:ROI_MethodAVG-F" = "ROI_Method")) %>%
  mutate(rule = recode(Interaction,
                       "BandNB:InverseLCMV" = "Band + Mask ~ ROI_Method",
                       "InverseLCMV:MaskTask-based" = "Mask + Band ~ ROI_Method",
                       "MaskTask-based:ROI_Method3SVD" = "Band + Inverse ~ ROI_Method",
                       "MaskTask-based:ROI_MethodAVG-F" = "Band + Inverse ~ ROI_Method")) %>%
  filter(`p-value` < p.threshold & value_col != '')
         

plotAndSaveInteraction <- function(value_col, group_col, other_col, rule) {
  if (group_col == "Inverse") {
    g1 <- "eLORETA"
    g2 <- "LCMV"
  } else if (group_col == "Mask") {
    g1 <- "Anatomical"
    g2 <- "Task-based"
  } else {
    error('Bad group column')
  }
  
  p_int <- plotProcessingInteraction(results, value_col, 
                                     group_col, g1, g2,
                                     as.formula(rule))
  
  savefile <- paste('interaction', value_col, group_col, 
                    other_col, sep='_')
  message(paste('plotAndSaveInteraction:', savefile))
  save_plot(file.path(output.folder, paste0(savefile, '.png')), 
            p_int, bg = "white", base_width = 8, base_height = 8.6)
  
  p_int
}

all_plots <- with(stats_int_plot, 
                  mapply(plotAndSaveInteraction, value_col, group_col,
                         other_col, rule, SIMPLIFY = F))

# Save a subset of plots to be used in the paper
p_snr_band_inverse <- all_plots[[1]] +
  xlab('SNR - eLORETA (dB)') + ylab('SNR - LCMV (dB)')
save_plot(file.path(output.folder, 'fig9supp1-snr-band-inverse-interaction.pdf'), 
          p_snr_band_inverse, base_width = 8, base_height = 8.6)
save_plot(file.path(output.folder, 'fig9supp1-snr-band-inverse-interaction.png'), 
          p_snr_band_inverse, bg = "white", base_width = 8, base_height = 8.6)

p_coh_within_mask_roi_method <- all_plots[[7]] +
  xlab('Within-hemisphere Coherence - Anatomical Definitions of ROIs') + 
  ylab('Within-hemisphere Coherence - Task-based Definitions of ROIs')
save_plot(file.path(output.folder, 'fig9supp2-coh-within-mask-roi-method-interaction.pdf'), 
          p_coh_within_mask_roi_method, base_width = 8, base_height = 8.6)
save_plot(file.path(output.folder, 'fig9supp2-coh-within-mask-roi-method-interaction.png'), 
          p_coh_within_mask_roi_method, bg = "white", base_width = 8, base_height = 8)

p_snr_inverse_band <- plotProcessingInteraction(results, 'LogSNR', 
                                                'Band', 'BB', 'NB',
                                                Inverse + Mask ~ ROI_Method) +
  xlab('SNR - Broadband (dB)') + 
  ylab('SNR - Narrowband (dB)')
save_plot(file.path(output.folder, 'fig9supp3-snr-inverse-band-interaction.pdf'), 
          p_snr_inverse_band, base_width = 8, base_height = 8.6)
save_plot(file.path(output.folder, 'fig9supp3-snr-inverse-band-interaction.png'), 
          p_snr_inverse_band, bg = "white", base_width = 8, base_height = 8)

### Save all the results
save(SNR_multiverse, pipelineEffects,
     tvals, pvals, file = file.path(r.path, 'pipeline_effects.RData'))