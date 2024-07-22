###
# Main script for the statistical analysis
###

## Set up the environment
library(corrplot)        # corrplot
library(data.table)      # rbindlist
library(dplyr)           # arrange
library(eegUtils)        # import.set, geom_topo
library(formula.tools)   # convert formula to string
library(lme4)            # lmer
library(lmerTest)        # Satterthwaite's method
library(kableExtra)      # cell_spec
library(magick)          # draw_image
library(ppcor)           # partial correlation
library(purrr)           # map
library(R.matlab)        # readMat
library(rmcorr)          # rmcorr
library(rstatix)         # cohens_d
library(testit)          # assert
library(tidyr)           # pivot_longer
library(scales)          # muted
library(stringr)         # str_replace
library(zeallot)         # %<-%
# Plotting
library(ggplot2)         # ggplot
library(ggh4x)           # multi-level facet grid plots (facet_nested)
library(ggsignif)        # significance annotations
library(ggrain)          # raincloud plot (geom_rain)
library(cowplot)         # ggdraw, aligning several ggplots

# Print the path to the working directory
(getwd())

# Path to the results exported from Matlab
deriv.path <- normalizePath('../data/derivatives/')
data.path <- file.path(deriv.path, 'task1')
eeg.path <- normalizePath('../data/preproc/task1/')
tex.path <- normalizePath('../results/tex/')
asset.path <- normalizePath('../assets/')

# Create an output folder for all results
r.path <- file.path(deriv.path, 'r')
dir.create(r.path)

# Create an output folder for all plots
plot.folder <- '../results/stats/'
dir.create(plot.folder)
plot.path <- normalizePath(plot.folder)

# Significance threshold
p.threshold <- 0.05
n.mc <- 6
p.mc <- p.threshold / n.mc

# Whether to run the code prepared during revision
revision <- FALSE

# Import the common functions
source('base.R')
source('plots.R')

# Define and export general settings (p-value thresholds and all formulas)
output_filename <- file.path(tex.path, 'general.tex')
source('formulas.R')

# Load the data for sanity checks and examples
sanity_filename <- file.path(data.path, 'BCI_MI_sanity_checks.mat')
c(sanity, conn_df) %<-% load_sanity(sanity_filename)

# Run analysis of performance
input_filename <- file.path(data.path, 'BCI_MI_performance_metadata_long.csv')
output_filename <- file.path(tex.path, 'metadata.tex')
source('metadata.R')

# Comparison of online and offline (CSP + rLDA) accuracy
input_filename <- file.path(data.path, 'BCI_MI_CSP_accuracy.mat')
output_filename <- file.path(tex.path, 'CSP_accuracy_AUC.tex')
source('CSP_accuracy_AUC.R')

# Make the pipeline overview plot
color_filename <- file.path(data.path, 'colors.mat')
source('pipeline_overview.R')

# Run analysis for Laplace SNR
input_filename <- file.path(data.path, 'BCI_MI_Laplace_SNR_results_long.mat')
output_filename <- file.path(tex.path, 'laplace_SNR.tex')
source('laplace_SNR.R')

# Run analysis for Laplace connectivity
input_filename <- file.path(data.path, 'BCI_MI_Laplace_connectivity_results_long.mat')
output_filename <- file.path(tex.path, 'laplace_connectivity.tex')
source('laplace_connectivity.R')

# Run analysis of multiverse results (resting-state)
input_filename <- file.path(data.path, 'BCI_MI_multiverse_rest_results_long.mat')
output_filename <- file.path(tex.path, 'multiverse_rest.tex')
source('multiverse_SNR.R')
source('multiverse_SNR_connectivity.R')
source('multiverse_connectivity_performance.R')
source('multiverse_connectivity_longitude.R')
source('multiverse_joint.R')

# Run analysis of pipeline-related effects
input_filename <- file.path(data.path, 'BCI_MI_multiverse_rest_results_long.mat')
output_filename <- file.path(tex.path, 'pipeline_effects.tex')
source('pipeline_effects.R')

# Revision
if (revision) {
  # Comparison of Fourier-based and Hilbert-based approaches for estimation of PS
  # revealed systematic differences that were confounding the effect of filtering
  source('fourier_vs_hilbert.R')
}
