%%% 
% Estimating the average SNR of mu oscillations at C3 and C4 electrodes
% after applying the Laplacian transform, plot FOOOF fits, export to R
%
% Output:
% 
% 1. BCI_MI_C3_C4_Laplace.mat - results_lap structure with fields:
%     analyzed_sessions [subjects x sessions] - 1 if the session was
%         analyzed, 0 otherwise
%     spec_laplace [subjects x sessions x laplace x freqs] - spectra of C3-
%         and C4-Laplace for each subject and session
%     freqs [freqs x 1] - frequencies corresponding to bins
%     snr_laplace [subjects x sessions x laplace] - SNR estimates for each
%         subject, session, and channel, NaN if not analyzed
%     r2_fooof [subjects x sessions x laplace] - r2 of the FOOOF fit for
%         each subject, session, and channel, NaN if not analyzed
%     fooof_fits {subjects x sessions x laplace} - fitted FOOOF models for
%         each subject, session, and channel, empty if not analyzed
%
% 2. BCI_MI_C3_C4_Laplace_results_long.mat - tabular data in the long
%    format for analysis in R (subject ID, session ID, channel (1: C3-Lap,
%    2: C4-Lap), accuracy, estimated SNR)
%%%

analyzed_sessions = zeros(n_subjects, n_sessions);
spec_laplace = zeros(n_subjects, n_sessions, n_laplace, n_freqs);
snr_laplace = NaN(n_subjects, n_sessions, n_laplace);
r2_fooof = NaN(n_subjects, n_sessions, n_laplace);
fooof_fits = cell(n_subjects, n_sessions, n_laplace);

%% Main loop
parfor (subject = 1:n_subjects, num_workers)
%     for subject = 1:n_subjects
    for session = 1:n_sessions
%         for session = 1
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        if ~exist([datapath filename], 'file')
            warning("Skipping subject %d session %d - file not found\n", subject, session);
            continue
        end

        % Load EEG
        EEG = pop_loadset('filepath', datapath, 'filename', filename);
        [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);

        % Extract data for condition - always use BB
        EEG_bb = extract_condition_EEG(EEG, EEG_narrow, {classes_to_analyze, 'rest', 'bb'}, events, windows);
        data_cnt = reshape(double(EEG_bb.data), n_sensors, []);

        % Apply Laplace transform
        data_lap_bb = zeros(n_laplace, size(data_cnt, 2));
        for ch = 1:n_laplace
            ind = laplace{ch, 1};
            ind_neighbors = laplace{ch, 2};
            data_lap_bb(ch, :) = data_cnt(ind, :) - mean(data_cnt(ind_neighbors, :), 1);
        end

        % Calculate the spectra
        spec_lap = pwelch(data_lap_bb', nfft, noverlap, nfft, EEG.srate)';
        spec_laplace(subject, session, :, :) = spec_lap;

        % Calculate the SNR
        for ch = 1:n_laplace
            % Run FOOOF
            fooof_results = fooof(freqs, squeeze(spec_lap(ch, :))', ...
                fit_range, fooof_settings, return_model);

            r2_fooof(subject, session, ch) = fooof_results.r_squared;
            snr_laplace(subject, session, ch) = ...
                calculate_SNR_fooof(fooof_results, methods, mu_band);
            fooof_fits{subject, session, ch} = fooof_results;
        end

        analyzed_sessions(subject, session) = 1;
    end
end

disp(sum(analyzed_sessions, [1 2]));

%% Save all results
results_lap = [];
results_lap.analyzed_sessions = analyzed_sessions;

% spectra
results_lap.spec_laplace = spec_laplace;
results_lap.freqs = freqs;
% SNR
results_lap.snr_laplace = snr_laplace;
results_lap.r2_fooof = r2_fooof;
results_lap.fooof_fits = fooof_fits;

save([savedata 'BCI_MI_C3_C4_Laplace.mat'], 'results_lap', '-v7.3');

%% Plot FOOOF fits
for subject = 1:n_subjects
    for session = 1:n_sessions
        if ~results_lap.analyzed_sessions(subject, session)
            continue;
        end

        h = figure('Position', [1 1 947 523]);
        for ch = 1:n_laplace
            chanlabel = all_chanlocs(laplace{ch, 1}).labels;
            ax = subplot(1, n_laplace, ch);
        
            fooof_results = fooof_fits{subject, session, ch};
            snr_mu = snr_laplace(subject, session, ch);
            display_SNR_fit_fooof(fooof_results, snr_mu, [], [], mu_band, ax);

            sgtitle(['S' num2str(subject) '/' num2str(session)]);
        end
        exportgraphics(h, [cfg.results.fooof_lap 'S' num2str(subject, '%02d') '_Session_' num2str(session, '%02d') '.png']);
        if (run_cmd)
            close(h);
        end

        fprintf('.');
    end
    fprintf('\n');
end

%% Export to R in the long format
load([savedata 'BCI_MI_task1_accuracy.mat']);   % task1_accuracy

% Export to R for statistical analysis
control_reshape = cell(n_subjects, n_sessions, n_laplace);
for subject = 1:n_subjects
    for session = 1:n_sessions
        for ch = 1:n_laplace
            if results_lap.analyzed_sessions(subject, session) && good_sessions(subject, session)
                control_reshape{subject, session, ch} = ...
                    [subject session ch task1_accuracy(subject, session)];
            else
                % mark session with NaN accuracy to be removed later
                control_reshape{subject, session, ch} = ...
                    [subject session ch NaN];
            end
        end
    end
end
metadata = cat(1, control_reshape{:});
snrdata = results_lap.snr_laplace;
data = cat(2, metadata, snrdata(:));
data(isnan(data(:, 4)), :) = [];
labels = {'Subject', 'Session', 'Channel', 'Accuracy', 'SNR'};
channels = {'C3-Lap', 'C4-Lap'};

save([savedata 'BCI_MI_C3_C4_Laplace_results_long.mat'], ...
    'data', 'labels', 'channels');