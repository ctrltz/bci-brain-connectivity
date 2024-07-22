###
# Joint multiverse analysis
#
# Outputs:
# 1. fig8-multiverse-joint-analysis.png
# 2. multiverse_rest.tex - table with statistical information (table S2)
# 3. multiverse_joint.RData - intermediate results of the analysis
###

# Prefix
prefix = "multiverseRest"

# Create output folder for images
output.folder <- file.path(plot.path, 'multiverse_joint')
dir.create(output.folder)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)

### Compile a summary of all the results of the multiverse analysis

all_results = rbindlist(c(
    list(
      SNR_vs_accuracy_within_results,
      SNR_vs_session_results
    ),
    SNR_vs_conn_results,
    conn_within_vs_acc_results,
    conn_across_vs_acc_results,
    conn_within_vs_session_results,
    conn_across_vs_session_results
  ),
  fill = T
) %>% as.data.frame()
names(all_results) <- c('Predictor', 'Response', 'Level', 
                        'SplitSignif', 'SplitSignifSameDir', 'Total', 
                        '<t-value>', '$\\beta$', 'df', 't-value', 'p-value', 
                        'CIMin', 'CIMax', 'Significant', 
                        'Significant.MC', 'Type', 'SNR Corr.')
all_results <- all_results %>% 
  mutate(Percent_Consistent = if_else(
    .data$Significant.MC,
    SplitSignifSameDir / Total,    # if significant, only count pipelines with the
                                   # same direction of the effect
    1 - SplitSignif / Total        # if non-significant, ignore the direction of the effect
  ))
all_results$CIDisp <- paste(
  "$[", format(round(all_results$CIMin, 3), nsmall = 3), ",", 
  format(round(all_results$CIMax, 3), nsmall = 3), "]$", sep = "")
all_results_sel <- all_results[, c('Predictor', 'Response',  
                                   '$\\beta$', 't-value', 'p-value', 
                                   'CIDisp', 'Percent_Consistent')]
sig.mc <- all_results %>%
  mutate_at(vars(Significant.MC), ~ if_else(., '*', '', missing = ''))

### Export it to LaTeX
all_results_latex <- all_results_sel %>%
  rename("Consistency" = "Percent_Consistent",
         "95\\% CI" = "CIDisp") %>%
  mutate_at(vars('Consistency'), round, 2) %>%
  mutate_at(vars('t-value', '$\\beta$', 'p-value'), round, 3) %>%
  mutate_at(vars('Predictor'), ~ case_when(
    . == 'LogSNR' ~ 'SNR',
    . == 'ImCoh_Within' ~ 'WH ImCoh',
    . == 'LagCoh_Within' ~ 'WH LagCoh',
    . == 'Coh_Within' ~ 'WH Coherence',
    . == 'ImCoh_Across' ~ 'AH ImCoh',
    . == 'LagCoh_Across' ~ 'AH LagCoh',
    . == 'Coh_Across' ~ 'AH Coherence',
    TRUE ~ .
  ))
# Format the p-values: bold denotes significance BEFORE correction for multiple
# comparison, stars - AFTER the correction
all_results_latex$`p-value` <- cell_spec(
  paste(if_else(all_results_latex[,'p-value'] < 0.001, '$<$0.001', 
                format(all_results_latex[,'p-value'], digits = 3)), 
        sig.mc[,'Significant.MC'], sep = ''), 
  format = 'latex', bold = all_results[, 'Significant'], escape = F)

all_results_latex %>%
  knitr::kable(format = "latex", booktabs = T, escape = F, align = 'c',
               linesep = c(
                 "\\midrule", "\\midrule",
                 "", "", "\\addlinespace", "", "", "\\midrule",
                 "", "", "\\addlinespace", "", "", "\\addlinespace", 
                 "", "", "\\addlinespace", "", "", "\\midrule",
                 "", "", "\\addlinespace", "", "", "\\addlinespace", 
                 "", "", "\\addlinespace", "", ""
               )) %>%
  tex.save(output_filename, "Summary", ., prefix = prefix)


### Plot an overview of the results of the joint analysis

new_order <- c(3:8, 1, 9:11, 15:17, 1, 12:14, 18:20, 
               2, 21:23, 27:29, 2, 24:26, 30:32)
reorder_pvals <- c(NA, all_results[new_order, 'p-value'])
reorder_pvals[c(15, 29)] <- NA

reorder_pc <- c(NA, all_results$Percent_Consistent[new_order])
reorder_pc[c(15, 29)] <- NA

joint_df <- data.frame(
  Category = rep(c('Average', rep('Within-hemisphere', times = 3),
                   rep('Across-hemisphere', times = 3)), times = 5),
  Measure  = rep(c('SNR', rep(c('ImCoh', 'LagCoh', 'Coherence'), times = 2)), times = 5),
  DummyX = rep('X', times = 35),
  DummyY = rep('Y', times = 35),
  Question = c(
    rep('Is Predicted by SNR', times = 7),
    rep('Predicts Accuracy', times = 7),
    rep('Predicts Accuracy\nafter Correction for SNR', times = 7),
    rep('Changes Longitudinally', times = 7),
    rep('Changes Longitudinally\nafter Correction for SNR', times = 7)
  ),
  Consistency = reorder_pc,
  p.value = reorder_pvals
)
joint_df <- joint_df[-c(1, 15, 29),]
joint_df$Significant.MC <- joint_df$p.value < p.mc
joint_df$DummyX <- factor(joint_df$DummyX)
joint_df$DummyY <- factor(joint_df$DummyY)
joint_df$Category <- factor(joint_df$Category, levels = c('Average', 'Within-hemisphere', 'Across-hemisphere'))
joint_df$Measure <- factor(joint_df$Measure, levels = c('SNR', 'ImCoh', 'LagCoh', 'Coherence'))
joint_df$Question <- factor(joint_df$Question, levels = c('Is Predicted by SNR', 'Predicts Accuracy', 
                                              'Predicts Accuracy\nafter Correction for SNR',
                                              'Changes Longitudinally',
                                              'Changes Longitudinally\nafter Correction for SNR'))


fig8 <- plotMultiverseJoint(joint_df, 'Consistency', 'Significant.MC', lim = 1, 
                            facet_rule = 'Question ~ Category + Measure')
ggsave(file.path(output.folder, 'fig8-multiverse-joint-analysis.png'), 
       fig8, bg = "white", width = 7, height = 5)


### Save all the results
save(all_results, joint_df,
     file = file.path(r.path, 'multiverse_joint.RData'))