%%%
% Configuration parameters for all the scripts
%%%

%% Task
task = 1;            % Task to be analyzed
                     % 1 - horizontal cursor control (left/right MI)
                     % 2 - vertical cursor control (rest vs both hands MI)
                     % 3 - 2D cursor control

%% Paths
datapath = [cfg.preproc.base 'task' num2str(task) '/'];                          % where is the preprocessed data located
savedata = [cfg.derivatives.base 'task' num2str(task) '/'];                      % where to store intermediate data
savepath = [cfg.results.base 'task' num2str(task) '/'];                          % where to store figures
savepath_csp = [savepath 'csp_patterns_sensor_space/'];                          % subfolder for figures with CSP patterns
savepath_group_csp = [savepath 'group_csp_patterns/'];                           % subfolder for figures with group CSP results
electrode_file = 'precomputed/standard-10-5-cap385_added_mastoids.elp';

mkdir_if_not_exists(datapath);
mkdir_if_not_exists(savedata);
mkdir_if_not_exists(savepath);
mkdir_if_not_exists(savepath_csp);
mkdir_if_not_exists(savepath_group_csp);

fprintf('Data path: %s\n\n', datapath);
fprintf('Derivatives path: %s\n\n', savedata);
fprintf('Figures path: %s\n\n', savepath);

%% Constants
n_subjects = 62;        % Total number of subjects in the dataset
n_sessions = 11;        % Max number of session per subject    

run_cmd = 1;            % Close figures automatically in command line mode
num_workers = 12;       % Number of workers for parallel computation

% Processing
downsample = 250;       % Downsample to 250 Hz
n_mirror_seconds = 8;   % Use mirroring to suppress edge effects due to filtering
n_mirror_pnts = n_mirror_seconds * downsample;   

% Events & windows of interest
events.target = '1';              % Target appears on the screen
events.fbend = '3';               % Feedback period ends
windows.rest = [-1.512 -0.012];   % Relative to the target onset
windows.prep = [0.492 1.992];     % Relative to the target onset
windows.fbend = [-1.512 -0.012];  % Relative to the trial end

% Frequency band to use for filtering and connectivity
mu_band = [9 15];       % Use the same band as for the online feedback

% PSD
nfft = 375;   % to match 1.5 s segments extracted from the data
srate = 250;  % data will be downsampled
noverlap = 0;
freqs = psdfreqvec('npts', nfft, 'Fs', srate, 'Range', 'half');
n_freqs = numel(freqs);
mu_band_bins = freqs >= mu_band(1) & freqs <= mu_band(2);

% CS
fres = floor(srate * 3 / 4);  % 0.67 Hz frequency resolution
band_freqbins = freqs >= mu_band(1) & freqs <= mu_band(2);
maxfreq = 30;                 % Store only up to 30 Hz
maxfreqbin = find(freqs <= maxfreq, 1, 'last');

% Laplace montage for C3, C4, CP3, and CP4
laplace = {
    26, [25, 27, 17, 35];    % 'C3',  {'C1', 'C5', 'FC3', 'CP3'};
    30, [29, 31, 21, 39];    % 'C4',  {'C2', 'C6', 'FC4', 'CP4'};
    35, [36, 34, 26, 44];    % 'CP3', {'CP1', 'CP5', 'C3', 'P3'};
    39, [38, 40, 30, 48];    % 'CP4', {'CP2', 'CP6', 'C4', 'P4'};
};
n_laplace = size(laplace, 1);

% FOOOF & SNR
methods = {'ratio'};
n_methods = numel(methods);
fit_range = [1 45];
fooof_freqbins = freqs >= fit_range(1) & freqs <= fit_range(2);
n_fooof_freqbins = sum(fooof_freqbins);
fooof_settings = struct('peak_width_limits', [2.0 12.0], 'max_n_peaks', 3);
return_model = 1;

% Examples of different SNR
subj_low_snr = 54;
subj_med_snr = 43;
subj_high_snr = 60;

% CSP
csp_make_plots = 1;
n_csp_comps = 6;
csp_source_threshold = 97.5;

% ROI
% Indices from the Harvard-Oxford atlas of:
%     'Left Precentral Gyrus'
%     'Left Postcentral Gyrus'
%     'Right Precentral Gyrus'
%     'Right Postcentral Gyrus'
ROI_inds = [7 17 55 65];
roi_labels = {'preL', 'postL', 'preR', 'postR'};
n_rois = numel(ROI_inds);

% Inverse
lambda = 0.05; % regularization parameter, 
               % larger means stronger regularization, smoother solution

% Classes
classes_to_analyze = [1, 2];
classes = {'right', 'left'};
n_classes = numel(classes_to_analyze);

% Time periods
periods = {'rest', 'prep', 'fbend'};
n_periods = numel(periods);

% Frequency bands
bands = {'bb', 'nb'};
n_bands = numel(bands);
