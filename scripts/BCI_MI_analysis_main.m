%%% 
% Main analysis script
%%%

%% Initialize the environment
% Path to the conda/python/any other environment with FOOOF installed
python_path = '/data/u_kapralov_software/miniconda3/envs/fooof-for-matlab/bin/python3';

% Setup all the paths for MATLAB scripts to work
cfg = BCI_MI_init_workspace({'tprod', 'haufe', 'eeglab', ...
    'bbci', 'fooof', 'python'}, python_path);

load(cfg.toolbox.btb);  % required for the BBCI toolbox
load cm17;              % colormaps for source space visualizations
rng('default');

%% Setup
BCI_MI_config;
% run_cmd = 0; num_workers = 0; % turn off parallel computation for debugging

%% Run all tests
% BCI_MI_run_tests();

%% Apply preprocessing using the saved information about bad data
% Load the file with all preprocessing info
load([cfg.preproc.aux 'BCI_MI_preproc_info.mat']);

% Apply preprocessing
BCI_MI_preprocessing;

%% Extract values of accuracy
BCI_MI_performance;

%% Prepare the eLORETA inverse operator (same for all subjects)
BCI_MI_eLORETA_signflip;
load([pwd() '/precomputed/BCI_MI_sa_eLoreta.mat']);

%% Prepare the LCMV inverse operator (separate for each session)
BCI_MI_LCMV;
load([pwd() '/precomputed/BCI_MI_LCMV.mat']);

%% Run group-level CSP analysis to identify task-relevant sources
BCI_MI_CSP_fit;

% Load results of group CSP and the voxel mask
load([savedata 'BCI_MI_group_CSP.mat']);
load([savedata 'BCI_MI_CSP_voxel_mask.mat']);
BCI_MI_CSP_apply;

% Load CSP spectra and prepare the figure
load([savedata 'BCI_MI_CSP_spectra.mat']);
BCI_MI_CSP_plot;

%% Run the analysis for C3/C4-Laplace
BCI_MI_C3_C4_Laplace;

%% Run the multiverse analysis
BCI_MI_multiverse_config;

% Calculate the spectra of ROI time series
BCI_MI_spectra_multiverse;

% Load the spectra to fit FOOOF and estimate the SNR
load([savedata 'BCI_MI_spectra_rois_multiverse.mat']);
BCI_MI_SNR_multiverse;

% Calculate the connectivity for all pipelines
BCI_MI_connectivity_multiverse;

% Export SNR and connectivity to R in the long format
load([savedata 'BCI_MI_connectivity_multiverse.mat']);   % coh & others
load([savedata 'BCI_MI_SNR_rois_multiverse.mat']);       % snr_*, r2_*
load([savedata 'BCI_MI_task1_accuracy.mat']);            % task1_accuracy
BCI_MI_multiverse_export;

%% Sanity check
BCI_MI_mu_power_contrasts;

BCI_MI_power_connectivity_spectra;

%% Export the parameters to LaTeX
BCI_MI_latex_export;

%% Make plots about the pipeline for the paper
load([savedata 'BCI_MI_spectra_rois_multiverse.mat']);
BCI_MI_pipeline_overview;
