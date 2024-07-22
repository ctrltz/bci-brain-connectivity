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
  levels(conn_df$Mask) <- c("Anatomical", "Task-based")

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
  levels(results$Mask) <- c("Anatomical", "Task-based")

  # Log-scale the SNR to make the distribution closer to normal and use dB
  results$LogSNR <- 10*log10(results$SNR)

  # Extract description of the pipelines
  pipeline_desc <- results %>%
    group_by(Pipeline) %>%
    summarise(
      across(c(Inverse, Band, ROI_Agg, ROI_Method, Mask),
             ~ .x[[1]]
      ))

  list(results, pipeline_desc)
}


### Statistical analyses


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


scaleAndFitLME <- function (x, formula, columns_to_scale, refit = F) {
  # Fit linear mixed model after scaling the columns
  #
  # Parameters:
  #   x - dataset
  #   formula - mixed model equation
  #   columns_to_scale - list of columns to scale
  #   refit - whether to try ML instead of REML if did not converge
  #
  # Returns:
  #   fm - fitted model
  
  # Scale values within each pipeline separately
  for (col in columns_to_scale) {
    x[col] <- scale(x[col])
  }

  # Fit the linear mixed model
  fm <- lmer(formula = formula, data = x, REML = T)
  if (refit && !has_converged(fm)) {
    message('LME did not converge, refitting without REML')
    fm <- lmer(formula = formula, data = x, REML = F)
  }
  fm
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


fitMultiverseSplitLME <- function(df, params, pipeline_desc, refit = F) {
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
            params$formula_split, params$cols.scale, refit = refit)
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


fitMultiverseJointLME <- function(df, params, refit = F) {
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
  fm <- scaleAndFitLME(df, params$formula_joint, params$cols.scale, refit = refit)
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
      split_significant = sum(split_sel$Significant.MC),
      # NOTE: the formula below is also used separately in two lines of multiverse_SNR.R
      split_significant_samedir = sum(split_sel$Significant.MC &
                                      (split_sel$t.value * joint_sel$t.value > 0)),
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


formula2tex <- function(formula, rename_mapping = c("LogSNR" = "SNR", 
                                                    "Band" = "Filt.",
                                                    "Inverse" = "Inv.",
                                                    "Mask" = "ROI", 
                                                    "ROI_Method" = "Extr.")) {
  # Export linear mixed model formulas to TeX while renaming some of the factors
  #
  # Parameters:
  #   formula - R formula object
  #   rename_mapping - named list with new factor names

  result <- as.character(formula)
  result <- gsub("~", "$\\sim$", result, fixed = T)
  result <- gsub("|", "$|$", result, fixed = T)
  result <- gsub("^2", "$^2$", result, fixed = T)
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


checkOutliers <- function(x, threshold = 4) {
  # Scale the data and check which elements are further than threshold SDs away
  # from the mean
  #
  # Parameters:
  #   x - value
  #   threshold - how many SDs it takes to be considered as an outlier
  #
  # Returns:
  #   T/F vector for each element (T = is an outlier)
  
  abs(scale(x) >= threshold)
}
