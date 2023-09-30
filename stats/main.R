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
library(scales)          # muted
library(stringr)         # str_replace
library(zeallot)         # %<-%
# Plotting
library(ggplot2)         # ggplot
library(ggh4x)           # multi-level facet grid plots (facet_nested)
library(ggrain)          # raincloud plot (geom_rain)
library(cowplot)         # ggdraw, aligning several ggplots

# Print the path to the working directory
(getwd())

# Path to the results exported from Matlab
data.path <- normalizePath('../data/preproc/replication/')
eeg.path <- normalizePath('../data/preproc/task1/')
tex.path <- normalizePath('../results/tex/')
asset.path <- normalizePath('../assets/')

# Create an output folder for all results
r.folder <- '../data/r/'
dir.create(r.folder)
r.path <- normalizePath('../data/r/')

# Create an output folder for all plots
plot.folder <- '../results/stats/'
dir.create(plot.folder)
plot.path <- normalizePath(plot.folder)

# Significance threshold
p.threshold <- 0.05
n.mc <- 6
p.mc <- p.threshold / n.mc

# Import the common functions
source('base.R')

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

# Make the pipeline overview plot
color_filename <- file.path(data.path, 'colors.mat')
source('pipeline_overview.R')

# Run analysis for C3/C4-Laplace
input_filename <- file.path(data.path, 'BCI_MI_C3_C4_Laplace_results_long.mat')
output_filename <- file.path(tex.path, 'C3_C4_Laplace.tex')
source('C3_C4_Laplace.R')

# Run analysis of multiverse results (resting-state)
input_filename <- file.path(data.path, 'BCI_MI_multiverse_rest_results_long.mat')
output_filename <- file.path(tex.path, 'multiverse_rest.tex')
prefix = "multiverseRest"
source('multiverse_SNR.R')
source('multiverse_connectivity_performance.R')
source('multiverse_connectivity_longitude.R')
source('multiverse_joint.R')

# Run analysis of pipeline-related effects
input_filename <- file.path(data.path, 'BCI_MI_multiverse_rest_results_long.mat')
output_filename <- file.path(tex.path, 'pipeline_effects.tex')
source('pipeline_effect.R')
