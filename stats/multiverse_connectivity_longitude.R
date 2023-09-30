###
# Longitudinal changes in connectivity
#
# Outputs:
# 1. fig6-multiverse-connectivity-performance-within.png
# 2. fig6supp1-multiverse-connectivity-performance-between-subject.png
# 3. multiverse_connectivity_longitude.RData - intermediate results of the 
#    analysis
###

### Longitudinal Changes in Connectivity Within Hemispheres

param_space <- list(
  list(formula_split = connectivity_vs_session_formula("ImCoh_Within"),
       formula_joint = connectivity_vs_session_formula("ImCoh_Within", joint = T),
       cols.scale = c('Session', 'ImCoh_Within'), 
       predictor = 'Session', response = 'WH ImCoh',
       measure.name = 'ImCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("LagCoh_Within"),
       formula_joint = connectivity_vs_session_formula("LagCoh_Within", joint = T),
       cols.scale = c('Session', 'LagCoh_Within'), 
       predictor = 'Session', response = 'WH LagCoh',
       measure.name = 'LagCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("Coh_Within"),
       formula_joint = connectivity_vs_session_formula("Coh_Within", joint = T),
       cols.scale = c('Session', 'Coh_Within'), 
       predictor = 'Session', response = 'WH Coherence',
       measure.name = 'Coherence', measure.type = 'Not Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("ImCoh_Within", SNRCorr = T),
       formula_joint = connectivity_vs_session_formula("ImCoh_Within", SNRCorr = T, joint = T),
       cols.scale = c('Session', 'LogSNR', 'ImCoh_Within'), 
       predictor = 'Session', response = 'WH ImCoh $|$ SNR',
       measure.name = 'ImCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("LagCoh_Within", SNRCorr = T),
       formula_joint = connectivity_vs_session_formula("LagCoh_Within", SNRCorr = T, joint = T),
       cols.scale = c('Session', 'LogSNR', 'LagCoh_Within'), 
       predictor = 'Session', response = 'WH LagCoh $|$ SNR',
       measure.name = 'LagCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("Coh_Within", SNRCorr = T),
       formula_joint = connectivity_vs_session_formula("Coh_Within", SNRCorr = T, joint = T),
       cols.scale = c('Session', 'LogSNR', 'Coh_Within'), 
       predictor = 'Session', response = 'WH Coherence $|$ SNR',
       measure.name = 'Coherence', measure.type = 'Corrected for SNR')
)

conn_within_vs_session_stats <- lapply(param_space, fitMultiverseSplitLME,
                                       df = results, pipeline_desc = pipeline_desc)
conn_within_vs_session_stats <- do.call(rbind, conn_within_vs_session_stats)
conn_within_vs_session_stats$Significant.MC <- conn_within_vs_session_stats$p.value < p.mc
assert("Within PS ~ Session [split] did not converge",
       all(conn_within_vs_session_stats$Converged))

conn_within_vs_session_joint_stats <- lapply(param_space, fitMultiverseJointLME,
                                             df = results)
conn_within_vs_session_joint_stats <- do.call(rbind, conn_within_vs_session_joint_stats)
assert("Within PS ~ Session [joint] did not converge",
       all(conn_within_vs_session_joint_stats$Converged))

conn_within_vs_session_results <- lapply(param_space, getSummary,
                                         split_stats = conn_within_vs_session_stats,
                                         joint_stats = conn_within_vs_session_joint_stats,
                                         commonFields = list(type2 = 'Within', 
                                                             level = 'Within'))
conn_within_vs_session_results <- lapply(conn_within_vs_session_results,
                                         \(x) {
                                           x$SNRCorr <- x$type
                                           x$type <- x$type2
                                           x$type2 <- NULL
                                           x
                                         })


### Longitudinal Changes in Connectivity Across Hemispheres
param_space <- list(
  list(formula_split = connectivity_vs_session_formula("ImCoh_Across"),
       formula_joint = connectivity_vs_session_formula("ImCoh_Across", joint = T),
       cols.scale = c('Session', 'ImCoh_Across'), 
       predictor = 'Session', response = 'AH ImCoh',
       measure.name = 'ImCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("LagCoh_Across"),
       formula_joint = connectivity_vs_session_formula("LagCoh_Across", joint = T),
       cols.scale = c('Session', 'LagCoh_Across'), 
       predictor = 'Session', response = 'AH LagCoh',
       measure.name = 'LagCoh', measure.type = 'Not Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("Coh_Across"),
       formula_joint = connectivity_vs_session_formula("Coh_Across", joint = T),
       cols.scale = c('Session', 'Coh_Across'), 
       predictor = 'Session', response = 'AH Coherence',
       measure.name = 'Coherence', measure.type = 'Not Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("ImCoh_Across", SNRCorr = T),
       formula_joint = connectivity_vs_session_formula("ImCoh_Across", SNRCorr = T, joint = T),
       cols.scale = c('Session', 'LogSNR', 'ImCoh_Across'), 
       predictor = 'Session', response = 'AH ImCoh $|$ SNR',
       measure.name = 'ImCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("LagCoh_Across", SNRCorr = T),
       formula_joint = connectivity_vs_session_formula("LagCoh_Across", SNRCorr = T, joint = T),
       cols.scale = c('Session', 'LogSNR', 'LagCoh_Across'), 
       predictor = 'Session', response = 'AH LagCoh $|$ SNR',
       measure.name = 'LagCoh', measure.type = 'Corrected for SNR'),
  list(formula_split = connectivity_vs_session_formula("Coh_Across", SNRCorr = T),
       formula_joint = connectivity_vs_session_formula("Coh_Across", SNRCorr = T, joint = T),
       cols.scale = c('Session', 'LogSNR', 'Coh_Across'), 
       predictor = 'Session', response = 'AH Coherence $|$ SNR',
       measure.name = 'Coherence', measure.type = 'Corrected for SNR')
)

conn_across_vs_session_stats <- lapply(param_space, fitMultiverseSplitLME,
                                       df = results, pipeline_desc = pipeline_desc)
conn_across_vs_session_stats <- do.call(rbind, conn_across_vs_session_stats)
conn_across_vs_session_stats$Significant.MC <- conn_across_vs_session_stats$p.value < p.mc
assert("Across PS ~ Session [split] did not converge",
       all(conn_across_vs_session_stats$Converged))

conn_across_vs_session_joint_stats <- lapply(param_space, fitMultiverseJointLME,
                                             df = results)
conn_across_vs_session_joint_stats <- do.call(rbind, conn_across_vs_session_joint_stats)
assert("Across PS ~ Session [joint] did not converge",
       all(conn_across_vs_session_joint_stats$Converged))

conn_across_vs_session_results <- lapply(param_space, getSummary,
                                         split_stats = conn_across_vs_session_stats,
                                         joint_stats = conn_across_vs_session_joint_stats,
                                         commonFields = list(type2 = 'Within', 
                                                             level = 'Within'))
conn_across_vs_session_results <- lapply(conn_across_vs_session_results,
                                         \(x) {
                                           x$SNRCorr <- x$type
                                           x$type <- x$type2
                                           x$type2 <- NULL
                                           x
                                         })


### Plot the Multiverse
lim_conn_session = ceilingn(max(abs(c(conn_within_vs_session_stats$Estimate,
                                      conn_across_vs_session_stats$Estimate))), 2)

conn_within_vs_session_stats <- within(conn_within_vs_session_stats, {
  fMeasure = factor(Measure, levels = c('ImCoh', 'LagCoh', 'Coherence'))
  Band = factor(Band, levels = c('BB', 'NB'), labels = c('Broadband', 'Narrowband'))
  fType = factor(Type, levels = c('Not Corrected for SNR', 'Corrected for SNR'))
})

p_conn_session_within <- plotMultiverseSplit(conn_within_vs_session_stats, 'Estimate', 'Significant.MC',
                                             lim = lim_conn_session, facet_rule = 'fType + Band ~ fMeasure + Mask')

conn_across_vs_session_stats <- within(conn_across_vs_session_stats, {
  fMeasure = factor(Measure, levels = c('ImCoh', 'LagCoh', 'Coherence'))
  Band = factor(Band, levels = c('BB', 'NB'), labels = c('Broadband', 'Narrowband'))
  fType = factor(Type, levels = c('Not Corrected for SNR', 'Corrected for SNR'))
})

p_conn_session_across <- plotMultiverseSplit(conn_across_vs_session_stats, 'Estimate', 'Significant.MC',
                                             lim = lim_conn_session, facet_rule = 'fType + Band ~ fMeasure + Mask')

fig6supp2 <- plot_grid(
  p_conn_session_within + theme(legend.position = 'none'),
  p_conn_session_across + theme(legend.position = 'none'),
  get_legend(p_conn_session_within),
  ncol = 1, rel_heights = c(1, 1, 0.2),
  labels = c('A', 'B', '')
)

save_plot(file.path(plot.path, prefix, 'fig6supp2-multiverse-connectivity-longitude.png'),
          fig6supp2, bg = "white", base_width = 8, base_height = 12)


### Save all the results
save(conn_within_vs_session_stats, conn_within_vs_session_joint_stats,
     conn_across_vs_session_stats, conn_across_vs_session_joint_stats,
     conn_within_vs_session_results, conn_across_vs_session_results,
     file = file.path(r.path, 'multiverse_connectivity_longitude.RData'))