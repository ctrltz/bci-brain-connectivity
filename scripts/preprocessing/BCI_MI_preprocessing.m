%%% 
% Apply preprocessing to the BCI MI dataset from (Stieger et al., 2021)
%
% Outputs:
% 1. results/tex/preprocessing.tex - list of channels available in the
% dataset
%%%

%% Load data & configure
% Structure par contains parameters to pass into the preprocessing function
par = cfg;
par.savepath = [cfg.results.base 'preprocessing/'];
par.savedata_raw_eeglab = [cfg.data.base 'eeglab/'];
par.electrode_file = electrode_file;

if ~exist(par.savepath, 'dir')
   mkdir(par.savepath);
   mkdir([par.savepath 'apply/']);
end

if ~exist(par.savedata_raw_eeglab, 'dir')
   mkdir(par.savedata_raw_eeglab);
end

par.run_cmd = run_cmd;
par.save_raw = 1;
par.save_preproc = 1;

%% Apply preprocessing to the original data
% Loop over the participants and apply preprocessing
parfor (subject = 1:n_subjects, num_workers)
% for subject = 1:n_subjects
    for session = 1:n_sessions
        prepinfo = preproc_info{subject, session};
        if ~isempty(prepinfo)
            BCI_MI_preprocessing_apply(subject, session, prepinfo, par);
        end
    end
end

%% Get the list of channels from an arbitrary session
% Original data is used - all channels are still there
subject = 1;
session = 1;
filename = ['S' num2str(subject) '_Session_' num2str(session) '.mat'];
ses_data = load([cfg.data.raw filename]);
chan_labels = ses_data.BCI.chaninfo.label;

% Fix the letter case for display purposes
chan_labels = cellfun(@(x) replace(replace(x, 'FP', 'Fp'), 'Z', 'z'), ...
    chan_labels, 'UniformOutput', false);

%% Export to LaTeX
fileID = fopen([cfg.results.tex 'preprocessing.tex'], 'w');
fprintf(fileID, '%% Data Format\n');
fprintf(fileID, '%s\n', export_to_latex("numChannelsOrig", numel(chan_labels), "%d"));
fprintf(fileID, '%s\n', export_to_latex("channelsOrig", strjoin(chan_labels, ", "), "%d"));
fclose(fileID);
