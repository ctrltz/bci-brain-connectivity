%%% 
% Preprocessing script for the BCI MI dataset from (Stieger et al., 2021)
% The script is provided for completeness and is not included in the main
% analysis since the results of preprocessing were validated and fixed
% manually if needed
%%%

%% Initialize environment
addpath('..');
cfg = BCI_MI_init_workspace({'eeglab', 'fastica'});
rng('default');

%% Load data & configure
% par accumulates optional preprocessing parameters, see description in the
% BCI_MI_preprocessing_fit.m
par = cfg;
par.savepath = [cfg.results.base 'preprocessing/'];
par.savedata_raw_eeglab = [cfg.data.base 'BCI_MI_EEGLAB/'];
par.electrode_file = electrode_file;

if ~exist(par.savepath, 'dir')
   mkdir(par.savepath)
   mkdir([par.savepath 'fit/']);
   mkdir([par.savepath 'apply/']);
end

par.run_cmd = 1;       % close figures automatically in command line mode
par.run_manual = 0;    % manually select bad trials, channels, and comps
par.save_raw = 0;      % save raw data in EEGLAB format
par.save_preproc = 1;  % save preprocessed data in EEGLAB format

% Preprocessing stage statuses
% (0 - stop, 1 - fit, 2 - apply, 3 - skip)
ss = [];
ss.stop = 0;
ss.fit = 1;
ss.apply = 2;
ss.skip = 3;
par.ss = ss;

% Preprocessing stages
stage.bad_channels_epochs = ss.fit;
stage.hp_averef = ss.apply;
stage.ica = ss.fit;
stage.ica_reject = ss.fit;
stage.peak2peak_reject = ss.skip;
par.stage = stage;

laplacian_channels = {'FC3', 'CP3', 'C1', 'C5', 'C3', 'FC4', 'CP4', 'C2', 'C6', 'C4'};
par.laplacian_channels = laplacian_channels;

num_workers = 12;

%% Fit preprocessing on the downsampled data
for subject = 1:n_subjects
    for session = 1:n_sessions
        BCI_MI_preprocessing_fit(subject, session, par);
    end
end

%% Edit preprocessing info manually to fix errors
subject = 1;
session = 1;
subjprefix = ['S' num2str(subject) '_Session_', num2str(session)];
preproc_info_filename = [subjprefix '_preproc_info.mat'];
tmp = load([cfg.preproc.aux preproc_info_filename]);
preproc_info = tmp.preproc_info;

% Add trials
trials_to_add = [1 2 3];
preproc_info.bad_trials = [preproc_info.bad_trials; trials_to_add']

% Add channels
channels_to_add = {'FP1', 'FP2'};
preproc_info.bad_channels = [preproc_info.bad_channels channels_to_add]

% Delete components
comps_to_reject = [1 2 3];
preproc_info.gcompreject(comps_to_reject) = 1;

% Save modified info
save([cfg.preproc.aux preproc_info_filename], 'preproc_info');

%% Merge all preprocessing info into one file
BCI_MI_preprocessing_merge;
