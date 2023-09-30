###
# Common functions used by other scripts
###


### Data loaders


load_sanity <- function (datapath) {
  # Load the sanity checks data (mu power contrasts, example spectra for SNR,
  # connectivity spectra)
  #
  # Parameters:
  #   datapath - path to file
  
  data <- readMat(datapath)
  conn_df <- as.data.frame(data$conn.data)
  colnames(conn_df) <- map(data$conn.labels, ~ .x[[1]])

  conn_df$Pipeline <- as.factor(conn_df$Pipeline)
  conn_df$Measure <- as.factor(conn_df$Measure)
  conn_df$Type <- as.factor(conn_df$Type)
  conn_df$Inverse <- as.factor(conn_df$Inverse)
  conn_df$BandViz <- 3 - conn_df$Band
  conn_df$Band <- as.factor(conn_df$Band)
  conn_df$Mask <- as.factor(conn_df$Mask)
  conn_df$ROI_Agg <- as.factor(conn_df$ROI_Agg)
  conn_df$ROI_Method <- as.factor(conn_df$ROI_Method)

  levels(conn_df$Pipeline) <- unlist(map(unlist(data$pipeline.labels), \(x) gsub('-normal', '', x, fixed = T)))
  levels(conn_df$Measure) <- c("ImCoh", "LagCoh", "Coherence")
  levels(conn_df$Type) <- c("Within-hemisphere", "Across-hemisphere")
  levels(conn_df$Inverse) <- c("eLORETA", "LCMV")
  levels(conn_df$Band) <- c("BB")
  levels(conn_df$ROI_Agg) <- c("1SVD", "3SVD", "AVG-F", "mask/1SVD", "mask/3SVD", "mask/AVG-F")
  levels(conn_df$ROI_Method) <- c("1SVD", "3SVD", "AVG-F")
  levels(conn_df$Mask) <- c("No Mask", "Mask")

  list(data, conn_df)
}


load_C3C4Laplace <- function (datapath) {
  # Load the SNR estimated at C3 and C4-Laplace
  #
  # Parameters:
  #   datapath - path to file
  
  data <- readMat(datapath)
  results <- as.data.frame(data$data)
  colnames(results) <- map(data$labels, ~ .x[[1]])

  results$Subject <- as.factor(results$Subject)
  results$Channel <- as.factor(results$Channel)
  levels(results$Channel) <- unlist(data$channels)

  # Log-scale the SNR to make the distribution closer to normal and use dB
  results$LogSNR <- 10*log10(results$SNR)
  results
}


load_multiverse <- function (datapath) {
  # Load SNR and connectivity estimated with different pipelines
  #
  # Parameters:
  #   datapath - path to file
  
  data <- readMat(datapath)
  results <- as.data.frame(data$data)
  colnames(results) <- map(data$labels, ~ .x[[1]])

  results$Pipeline <- as.factor(results$Pipeline)
  results$Subject <- as.factor(results$Subject)
  results$Inverse <- as.factor(results$Inverse)
  # make broadband value larger than narrowband to map to linewidth later
  results$BandViz <- 3 - results$Band
  results$Band <- as.factor(results$Band)
  results$Mask <- as.factor(results$Mask)
  results$ROI_Agg <- as.factor(results$ROI_Agg)
  results$ROI_Method <- as.factor(results$ROI_Method)

  levels(results$Pipeline) <- unlist(map(unlist(data$multiverse.labels), \(x) gsub('-normal', '', x, fixed = T)))
  levels(results$Inverse) <- c("eLORETA", "LCMV")
  levels(results$Band) <- c("BB", "NB")
  levels(results$ROI_Agg) <- c("1SVD", "3SVD", "AVG-F", "mask/1SVD", "mask/3SVD", "mask/AVG-F")
  levels(results$ROI_Method) <- c("1SVD", "3SVD", "AVG-F")
  levels(results$Mask) <- c("No Mask", "Mask")

  # Log-scale the SNR to make the distribution closer to normal and use dB
  results$LogSNR <- 10*log10(results$SNR)

  # Select only broad-band pipelines for analyses that include only SNR and no connectivity
  # since SNR was computed for broad-band data and copied to narrow-band analogues.
  results_SNR <- results[grepl("/bb/", results$Pipeline, fixed = T),]
  results_SNR$Pipeline <- factor(results_SNR$Pipeline)

  # Extract description of the pipelines
  pipeline_desc <- results %>%
    group_by(Pipeline) %>%
    summarise(
      across(c(Inverse, Band, ROI_Agg, ROI_Method, Mask),
             ~ .x[[1]]
      ))

  list(results, results_SNR, pipeline_desc)
}


### Statistical analyses


scaleAndFitLME <- function (x, formula, columns_to_scale,
                            REML = T) {
  # Fit linear mixed model after scaling the columns
  #
  # Parameters:
  #   x - dataset
  #   formula - mixed model equation
  #   columns_to_scale - list of columns to scale
  #   REML - whether to use REML criterion for optimization
  #
  # Returns:
  #   fm - fitted model
  
  # Scale values within each pipeline separately
  for (col in columns_to_scale) {
    x[col] <- scale(x[col])
  }

  # Fit the linear mixed model
  fm <- lmer(formula = formula, data = x, REML = REML)
  fm
}


has_converged <- function(fm) {
  # Helper function to check that LME converged
  #
  # Parameters:
  #   fm - fitted model
  #
  # Returns:
  #   TRUE if converged, FALSE if not
  #
  # Source: https://stackoverflow.com/a/72128391
  
  retval <- is.null(unlist(fm@optinfo$conv$lme4))
  if (retval) {
    retval
  } else {
    isSingular(fm)
  }
}

restoreBetas <- function(beta0, beta1, data, x, y) {
  # Convert betas calculated using standardized data back to original scale
  #
  # Parameters:
  #   beta0 - intercept
  #   beta1 - slope
  #   data - dataset
  #   x - x column before scaling
  #   y - y column before scaling
  #
  # Returns:
  #   intercept and slope converted to the original scale
  
  # Calculate mean and SD of the original data
  # NOTE: as.matrix is used to make it work with tibbles as well
  muX = mean(as.matrix(data[,x]))
  muY = mean(as.matrix(data[,y]))
  sigmaX = sd(as.matrix(data[,x]))
  sigmaY = sd(as.matrix(data[,y]))

  # Linear regression for standardized data:
  # (Y - muY) / sigmaY = beta0 + beta1 * (X - muX) / sigmaX
  # Y = muY + sigmaY * beta0 + beta1 * sigmaY * (X - muX) / sigmaX
  # Y = muY + sigmaY * (beta0 - beta1 * muX / sigmaX) + beta1 * sigmaY / sigmaX * X
  #     ---------------------------------------------   -----------------------
  #                    beta0_orig                             beta1_orig
  beta0_orig <- muY + sigmaY * (beta0 - beta1 * muX / sigmaX)
  beta1_orig <- beta1 * sigmaY / sigmaX

  c(beta0_orig, beta1_orig)
}


corrSubjectAverage <- function (df, x, y, group) {
  # Estimate between-subject effects (first average all values for the subject, 
  # then correlate)
  #
  # Parameters:
  #   df - data frame
  #   x, y - columns to correlate
  #   group - values are averaged within the group column before correlating
  #
  # Returns:
  #   correlation test between group-averaged df$x and df$y
  
  # Get subject-level means
  data <- aggregate(df[, c(x, y)], list(df[, group]), mean)

  # Run correlation test
  cor.test(data[,x], data[,y])
}


partialCorrSubjectAverage <- function (df, x, y, z, group) {
  # Estimate between-subject effects controlling for a third factor (first 
  # average all values for the subject, then use partial correlation)
  #
  # Parameters:
  #   df - data frame
  #   x, y - columns to correlate
  #   z - column to account for in the partial correlation
  #   group - values are averaged within the group column before correlating
  #
  # Returns:
  #   partial correlation test between group-averaged df$x and df$y controlling
  #   for group-averaged df$z
  
  # Get subject-level means
  data <- aggregate(df[, c(x, y, z)], list(df[, group]), mean)

  # Run correlation test
  pcor.test(data[,x], data[,y], data[,z])
}


fitMultiverseSplitBetween <- function(df, params, pipeline_desc) {
  # Estimate between-subject effects in the multiverse analysis
  #
  # Parameters:
  #   df - dataset
  #   params - list of parameters for corrSubjectAverage (x, y, group) as well 
  #     as some metadata for future joint analysis (measure name and type)
  #   pipeline_desc - dataframe with the description of all pipelines
  #
  # Returns:
  #   dataframe with the estimated correlations, t-values, degrees of freedom, 
  #     and p-values
  
  message(paste("fitMultiverseSplitBetween:", params$x, '~', params$y))
  cts <- by(df, df$Pipeline, corrSubjectAverage,
            params$x, params$y, params$group)
  coefs <- do.call(rbind, lapply(cts,
                                 \(ct) unlist(ct[c('statistic', 'parameter', 'p.value', 'estimate')]))) %>%
    as.data.frame() %>%
    rename("t.value" = "statistic.t", "df" = "parameter.df", "Estimate" = "estimate.cor") %>%
    tibble::rownames_to_column("Pipeline") %>%
    inner_join(., pipeline_desc, by = "Pipeline")
  coefs$Measure <- params$measure.name
  coefs$Type <- params$measure.type
  coefs$Significant <- coefs$p.value < p.threshold

  coefs
}


fitMultiverseSplitBetweenPartial <- function(df, params, pipeline_desc) {
  # Estimate between-subject effects in the multiverse analysis controlling for
  # a factor using partial correlation
  #
  # Parameters:
  #   df - dataset
  #   params - list of parameters for partialCorrSubjectAverage (x, y, z, group)
  #     as well as some metadata for future joint analysis (measure name / type)
  #   pipeline_desc - dataframe with the description of all pipelines
  #
  # Returns:
  #   dataframe with the estimated correlations, t-values, degrees of freedom, 
  #     and p-values
  
  message(paste('fitMultiverseSplitBetweenPartial:', 
                params$x, '~', params$y, '|', params$z))
  cts <- by(df, df$Pipeline, partialCorrSubjectAverage,
            params$x, params$y, params$z, params$group)
  coefs <- do.call(rbind, lapply(cts,
                                 \(ct) unlist(ct[c('statistic', 'n', 'p.value', 'estimate')]))) %>%
    as.data.frame() %>%
    rename("t.value" = "statistic", "df" = "n", "Estimate" = "estimate") %>%
    tibble::rownames_to_column("Pipeline") %>%
    inner_join(., pipeline_desc, by = "Pipeline")
  coefs$Measure <- params$measure.name
  coefs$Type <- params$measure.type
  coefs$Significant <- coefs$p.value < p.threshold

  coefs
}


fitMultiverseSplitLME <- function(df, params, pipeline_desc, REML = T) {
  # Estimate within-subject effects in the split multiverse analysis
  #
  # Parameters:
  #   df - dataset
  #   params - list of parameters for scaleAndFitLME (columns_to_scale, formula)
  #     and some metadata for future joint analysis (measure name and type)
  #   pipeline_desc - dataframe with the description of all pipelines
  #   REML - whether to use REML criterion for LME optimization
  #
  # Returns:
  #   dataframe with the estimated betas, t-values, degrees of freedom, 
  #     and p-values
  
  message(paste('fitMultiverseSplitLME:', params$formula_split))
  fms <- by(df, df$Pipeline, scaleAndFitLME,
            params$formula_split, params$cols.scale, REML = REML)
  coefs <- do.call(rbind, lapply(fms,
                                 \(fm) unlist(lmerTest:::get_coefmat(fm)[params$predictor,]))) %>%
    as.data.frame() %>%
    rename("p.value" = "Pr(>|t|)", "t.value" = "t value", "SE" = "Std. Error") %>%
    tibble::rownames_to_column("Pipeline") %>%
    inner_join(., pipeline_desc, by = "Pipeline")
  coefs$Measure <- params$measure.name
  coefs$Type <- params$measure.type
  coefs$Significant <- coefs$p.value < p.threshold
  coefs$Converged <- unlist(lapply(fms, has_converged))
  coefs$Singular <- unlist(lapply(fms, isSingular))

  coefs
}


fitMultiverseJointLME <- function(df, params) {
  # Estimate within-subject effects in the joint multiverse analysis
  #
  # Parameters:
  #   df - data frame
  #   params - parameters for scaleAndFitLME (formula, columns to scale) and
  #     metadata (measure name and type)
  #
  # Returns:
  #   dataframe with the estimated beta, t-values, degrees of freedom, CI,
  #     and p-values
  
  message(paste("fitMultiverseJointLME:", params$formula_joint))
  fm <- scaleAndFitLME(df, params$formula_joint, params$cols.scale)
  ci <- confint(fm)
  coefs <- unlist(lmerTest:::get_coefmat(fm)[params$predictor,]) %>%
    t() %>% as.data.frame() %>%
    rename("p.value" = "Pr(>|t|)", "t.value" = "t value", "SE" = "Std. Error")
  coefs$Measure <- params$measure.name
  coefs$Type <- params$measure.type
  coefs$CIMin <- ci[[params$predictor, 1]]
  coefs$CIMax <- ci[[params$predictor, 2]]
  coefs$Significant <- coefs$p.value < p.threshold
  coefs$Significant.MC <- coefs$p.value < p.mc
  coefs$Converged <- has_converged(fm)
  coefs$Singular <- isSingular(fm)

  coefs
}


getSummary <- function(params, split_stats, joint_stats, commonFields) {
  # Organize the results of split and joint multiverse analyses
  #
  # Parameters:
  #   params - params provided to the LME fit functions
  #   split_stats - aggregated results of the split analysis
  #   joint_stats - aggregated results of the joint analysis
  #   commonFields - fields to add to the summary
  #
  # Returns:
  #   summary of the results
  
  split_sel <- split_stats[split_stats$Measure == params$measure.name &
                           split_stats$Type == params$measure.type, ]
  joint_sel <- joint_stats[joint_stats$Measure == params$measure.name &
                           joint_stats$Type == params$measure.type, ]

  c(
    list(
      predictor = params$predictor,
      response = params$response,
      type = params$measure.type,
      split_significant = sum(split_sel$Significant),
      split_total = length(split_sel$Significant),
      split_t.value = mean(split_sel$t.value),
      joint_estimate = joint_sel$Estimate,
      joint_df = joint_sel$df,
      joint_t.value = joint_sel$t.value,
      joint_p.value = joint_sel$p.value,
      joint_ci_min = joint_sel$CIMin,
      joint_ci_max = joint_sel$CIMax,
      joint_result = joint_sel$p.value < p.threshold,
      joint_result.MC = joint_sel$p.value < p.mc
    ),
    commonFields
  )
}


###
# Plotting
###


plotPipelineHistograms <- function(data, param, nrow, plot.path, prefix, width, height) {
  # Density plots of estimated values for different pipelines
  #
  # Parameters:
  #   data - data frame
  #   param - list with variable to plot (param$x) and filename (param$name)
  #   nrow - number of rows for facet wrap
  #   plot.path - where to save the plots
  #   prefix - prefix of the filename
  #   width, height - parameters for ggsave
  
  p <- ggplot(data, mapping = aes(x = .data[[param$x]])) +
    geom_histogram(mapping = aes(y = after_stat(density)),
                   bins = 50, alpha = 0.5) +
    geom_density() +
    facet_wrap(. ~ Pipeline, nrow = nrow) +
    plot_theme

  output_filename <- file.path(plot.path,
                               paste(prefix, param$name, sep = ''))
  ggsave(output_filename, p, width = width, height = height)
}


plotMultiverseBase <- function(df, val, sig, rows, cols, facet_rule,
                               val.name = 'Estimate', sig.name = 'Significant', 
                               lim = NULL, shapes = c("TRUE" = 16, "FALSE" = 1),
                               fontface = c('bold', 'plain'), strip.y.angle = 90) {
  # Plot the results of the multiverse analysis as a table for facilitating
  # visual comparison of results of different pipelines. Color codes the 
  # continuous value of interest (e.g., beta, t or correlation coefficient), 
  # dots indicate the significance of the effect of interest.
  #
  # The idea behind this function is to use nested rows and columns for coding
  # different pipeline steps. It is expected that the pipeline has at least
  # two steps (x and y), additional steps might be provided as a facet_rule
  #
  # Example: 
  #   pipeline: A -> B -> C -> D
  #   desired output: rows code B1, ..., Bn, grouped by A1, ..., An
  #                   columns code D1, ..., Dn, grouped by C1, ..., Cn
  #   arguments: x = D, y = B, facet_rule = A ~ C
  #
  # Parameters:
  #   df - data frame
  #   val - column to use for displaying the value of interest
  #   sig - column to use for displaying significance of the results
  #   rows, cols - pipeline steps (lowest level of hierarchy)
  #   facet_rule - pipeline steps (higher levels of hierarchy, 
  #                                cols3 + cols2 ~ rows3 + rows2)
  #   val.name - colorbar title for the value
  #   sig.name - legend title for the significance
  #   lim - data limit (calculated automatically if not provided)
  #   shapes - how to display significance
  #   fontface - optional styling for nested pipeline levels
  #   strip.y.angle - optional rotation of y labels
  # 
  # Returns:
  #   ggplot
  
  # Calculate data limits if not provided manually
  if (is.null(lim)) {
    lim = ceilingn(max(abs(df[,val])), 2)
  }

  # Calculate the number of terms in both sides of the facet rule
  facet_formula <- as.formula(facet_rule)
  parts <- unlist(strsplit(as.character(facet_formula), ' ~ '))
  n_lhs <- length(unlist(strsplit(parts[1], '[+]')))
  n_rhs <- length(unlist(strsplit(parts[2], '[+]')))

  # Prepare the main plot, colormap can be adapted separately
  p <- ggplot(df, aes(x = .data[[cols]], y = .data[[rows]])) +
    geom_raster(aes(fill = .data[[val]])) +
    geom_point(aes(shape = .data[[sig]]), color = "black",
               size = 3, fill = NA, stroke = 1) +
    scale_shape_manual(labels = c('TRUE' = 'Yes', 'FALSE' = 'No'),
                       limits = c('TRUE', 'FALSE'), values = shapes, drop = FALSE,
                       guide = guide_legend(ncol = 1, order = 2,
                                            title = sig.name, title.position = 'top')) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    guides(fill = guide_colorbar(order = 1, title.position = 'top')) +
    facet_nested(facet_formula, switch = "y",
                 strip = strip_nested(
                   text_x = elem_list_text(face = tail(fontface, n_rhs)),
                   text_y = elem_list_text(face = tail(fontface, n_lhs),
                                           angle = strip.y.angle),
                   by_layer_x = T, by_layer_y = T
                 )) +
    labs(fill = val.name) +
    theme_classic() +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      strip.placement = "outside"
    )

  p
}


plotMultiverseSplit <- function(df, val, sig, facet_rule, 
                                val.name = 'Estimate', lim = NULL) {
  # Plot the results of the split multiverse analysis
  # (setting up the appearance)
  #
  # Parameters:
  #   df - data frame
  #   val - column to use for displaying the value of interest
  #   sig - column to use for displaying significance of the results
  #   facet_rule - pipeline steps (higher levels of hierarchy, 
  #                                cols3 + cols2 ~ rows3 + rows2)
  #   val.name - legend title for the value of interest
  #   lim - data limit (calculated automatically if not provided)
  #
  # Returns:
  #   ggplot
  
  if (is.null(lim)) {
    lim = ceilingn(max(abs(df[,val])), 2)
  }
  
  p <- plotMultiverseBase(df, val, sig, 'Inverse', 'ROI_Method', 
                          facet_rule, val.name = val.name, lim = lim)
  
  p + scale_fill_gradient2(low = muted("blue"), high = muted("red"),
                           limits = c(-lim, lim), breaks = c(-lim, 0, lim)) +
    guides(y = guide_axis(angle = 90))
}


plotMultiverseJoint <- function(df, val, sig, facet_rule, lim = NULL) {
  # Plot the results of the joint multiverse analysis
  # (setting up the appearance)
  #
  # Parameters:
  #   df - data frame
  #   val - column to use for displaying the value of interest
  #   sig - column to use for displaying significance of the results
  #   facet_rule - pipeline steps (higher levels of hierarchy, 
  #                                cols3 + cols2 ~ rows3 + rows2)
  #   lim - data limit (calculated automatically if not provided)
  #
  # Returns:
  #   ggplot
  
  p <- plotMultiverseBase(df, val, sig, 'DummyY', 'DummyX', facet_rule, 
                          lim = lim, val.name = 'Consistency', strip.y.angle = 0)

  p + scale_fill_gradient2(low = "#FD8060", mid = "#FEE191", high = "#B0D8A4", 
                           midpoint = 0.5, limits = c(0, lim)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
}


###
# Export to TeX
###


format.pval <- function(p.value, threshold, digits) {
  # Format p-value before exporting (leave as is or replace with "< threshold")
  #
  # Parameters:
  #   p.value - actual p-value
  #   threshold - threshold value (e.g., 0.001)
  #   digits - number of significant digits to keep in the output
  #
  # Returns:
  #   formatted p-value

  if (p.value < threshold) {
    paste("<", format(threshold, digits = digits))
  } else {
    paste("=", format(p.value, digits = digits))
  }
}


formula2tex <- function(formula, rename_mapping = c("LogSNR" = "SNR", "Inverse" = "Inv.",
                                                    "ROI_Method" = "ROI")) {
  # Export linear mixed model formulas to TeX while renaming some of the factors
  #
  # Parameters:
  #   formula - R formula object
  #   rename_mapping - named list with new factor names

  result <- as.character(formula)
  result <- gsub("~", "$\\sim$", result, fixed = T)
  result <- gsub("|", "$|$", result, fixed = T)
  result <- str_replace_all(result, rename_mapping)

  result
}


tex.save <- function (filename, varname, value = NULL, prefix = '') {
  # Export value as a newcommand macro for TeX
  # \newcommand{prefixvarname}{value}
  #
  # Parameters:
  #   filename - file to export to
  #   varname - name of the newcommand or a comment string
  #   value - value to export or NULL if varname is a comment
  #   prefix - prefix for the newcommand name
  
  f <- file(filename, "a")
  
  if (is.null(value)) {
    writeLines(varname, con = f)
  } else {
    writeLines(sprintf("\\newcommand{\\%s%s}{%s}", prefix, varname, value),
               con = f);
  }
  close(f);
}


###
# Utility functions
###

ceilingn <- function(x, digits) {
  # Round the number up to the nearest decimal that has a given number of 
  # floating point digits
  #
  # Parameters:
  #   x - value
  #   digits - number of digits to keep
  #
  # Returns:
  #   the rounded value
  
  ceiling(x * 10 ^ digits) / 10 ^ digits
}
