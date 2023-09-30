###
# Equations for all linear mixed models
#
# Outputs:
# 1. general.tex - equations and p-thresholds exported to TeX
###


# Prefix for TeX output
prefix = "formula"

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)

# Export actual threshold of significance for p-values
tex.save(output_filename, "\n% Significance Testing\n")
tex.save(output_filename, "pSignificant", 
         format(p.threshold), prefix = '')
tex.save(output_filename, "numComparisons",
         format(n.mc), prefix = '')

# Longitudinal changes in accuracy
accuracy_vs_session_formula <- Accuracy ~ 1 + Session + (Session | Subject)

tex.save(output_filename, "\n% Formulas\n")
tex.save(output_filename, "AccuracySession",
         formula2tex(accuracy_vs_session_formula), prefix = prefix)

# Effect of SNR on accuracy
snr_vs_accuracy_formula <- Accuracy ~ 1 + LogSNR + (LogSNR|Subject)
snr_vs_accuracy_formula_joint <- update.formula(snr_vs_accuracy_formula, 
                                                ~ . + (1|Pipeline))

tex.save(output_filename, "AccuracySNR",
         formula2tex(snr_vs_accuracy_formula), prefix = prefix)

# Longitudinal changes in SNR
session_vs_snr_formula <- LogSNR ~ 1 + Session + (Session | Subject)
session_vs_snr_formula_joint <- update.formula(session_vs_snr_formula, 
                                               ~ . + (1|Pipeline))

tex.save(output_filename, "SNRSession",
         formula2tex(session_vs_snr_formula), prefix = prefix)

# Effect of SNR on connectivity
snr_vs_connectivity_formula_base <- PS ~ 1 + LogSNR + (1|Subject)

snr_vs_connectivity_formula <- function(conn, joint = F) {
  # Dynamically construct PS ~ SNR formula for different PS measures
  #
  # Parameters:
  #   conn - column with the PS measure of interest
  #   joint - add (1|Pipeline) as a random factor if TRUE
  #
  # Returns:
  #   constructed formula
  
  f_base <- as.character(snr_vs_connectivity_formula_base)
  f <- as.formula(str_replace_all(f_base, c("PS" = conn)))

  if (joint) {
    f <- update.formula(f, ~ . + (1|Pipeline))
  }

  f
}

tex.save(output_filename, "SNRConnectivity",
         formula2tex(snr_vs_connectivity_formula_base), prefix = prefix)

# Effect of SNR on connectivity
accuracy_vs_connectivity_formula_base <- Accuracy ~ 1 + PS + (1|Subject)
accuracy_vs_connectivity_formula_SNRCorr <- Accuracy ~ 1 + LogSNR + PS + (1|Subject)

accuracy_vs_connectivity_formula <- function(conn, SNRCorr = F, joint = F) {
  # Dynamically construct Accuracy ~ PS formula for different PS measures
  #
  # Parameters:
  #   conn - column with the PS measure of interest
  #   SNRCorr - include SNR in the model as a covariate
  #   joint - add (1|Pipeline) as a random factor if TRUE
  #
  # Returns:
  #   constructed formula
  
  base_formula <- accuracy_vs_connectivity_formula_base
  if (SNRCorr) {
    base_formula <- accuracy_vs_connectivity_formula_SNRCorr
  }

  f_base <- as.character(base_formula)
  f <- as.formula(str_replace_all(f_base, c("PS" = conn)))
  
  if (joint) {
    f <- update.formula(f, ~ . + (1|Pipeline))
  }

  f
}

tex.save(output_filename, "AccuracyConnectivity",
         formula2tex(accuracy_vs_connectivity_formula_base), prefix = prefix)
tex.save(output_filename, "AccuracyConnectivitySNR",
         formula2tex(accuracy_vs_connectivity_formula_SNRCorr), prefix = prefix)

# Longitudinal changes in connectivity
connectivity_vs_session_formula_base <- PS ~ 1 + Session + (1|Subject)
connectivity_vs_session_formula_SNRCorr <- PS ~ 1 + LogSNR + Session + (1|Subject)

connectivity_vs_session_formula <- function(conn, SNRCorr = F, joint = F) {
  # Dynamically construct PS ~ Session formula for different PS measures
  #
  # Parameters:
  #   conn - column with the PS measure of interest
  #   SNRCorr - include SNR in the model as a covariate
  #   joint - add (1|Pipeline) as a random factor if TRUE
  #
  # Returns:
  #   constructed formula
  
  base_formula <- connectivity_vs_session_formula_base
  if (SNRCorr) {
    base_formula <- connectivity_vs_session_formula_SNRCorr
  }

  f_base <- as.character(base_formula)
  f <- as.formula(str_replace_all(f_base, c("PS" = conn)))

  if (joint) {
    f <- update.formula(f, ~ . + (1|Pipeline))
  }

  f
}

tex.save(output_filename, "ConnectivitySession",
         formula2tex(connectivity_vs_session_formula_base), prefix = prefix)
tex.save(output_filename, "ConnectivitySessionSNR",
         formula2tex(connectivity_vs_session_formula_SNRCorr), prefix = prefix)

# Effect of different processing methods on SNR
snr_vs_processing_formula <- LogSNR ~ Inverse + ROI_Method + Mask + (1|Subject) + (1|Pipeline)

tex.save(output_filename, "SNRProcessing",
         formula2tex(snr_vs_processing_formula), prefix = prefix)

# Effect of different processing methods on connectivity
connectivity_vs_processing_formula_base <- PS ~ Inverse + ROI_Method + Band + Mask + (1|Subject) + (1|Pipeline)
connectivity_vs_processing_formula_SNRCorr <- PS ~ LogSNR + Inverse + ROI_Method + Band + Mask + (1|Subject) + (1|Pipeline)

connectivity_vs_processing_formula <- function(conn, SNRCorr = F) {
  # Dynamically construct PS ~ Processing formula for different PS measures
  #
  # Parameters:
  #   conn - column with the PS measure of interest
  #   SNRCorr - include SNR in the model as a covariate
  #
  # Returns:
  #   constructed formula
  
  base_formula <- connectivity_vs_processing_formula_base
  if (SNRCorr) {
    base_formula <- connectivity_vs_processing_formula_SNRCorr
  }

  f_base <- as.character(base_formula)
  f <- as.formula(str_replace_all(f_base, c("PS" = conn)))

  f
}

tex.save(output_filename, "ConnectivityProcessing",
         formula2tex(connectivity_vs_processing_formula_base), prefix = prefix)
tex.save(output_filename, "ConnectivityProcessingSNR",
         formula2tex(connectivity_vs_processing_formula_SNRCorr), prefix = prefix)
