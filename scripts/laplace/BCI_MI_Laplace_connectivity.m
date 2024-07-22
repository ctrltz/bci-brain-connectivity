%%% 
% Estimating the within- and across-hemisphere PS using C3, C4, CP3, CP4
% after applying the Laplacian transform, export to R
%
% Output:
% 
% 1. BCI_MI_Laplace_connectivity.mat - results_lap_conn structure with fields:
%     analyzed_sessions [subjects x sessions] - 1 if the session was
%         analyzed, 0 otherwise
%     spec_laplace [subjects x sessions x laplace x freqs] - spectra of all
%         Laplace-filtered channels for each subject and session
%     freqs [freqs x 1] - frequencies corresponding to bins
%     snr_laplace [subjects x sessions x laplace] - SNR estimates for each
%         subject, session, and channel, NaN if not analyzed
%     r2_fooof [subjects x sessions x laplace] - r2 of the FOOOF fit for
%         each subject, session, and channel, NaN if not analyzed
%     fooof_fits {subjects x sessions x laplace} - fitted FOOOF models for
%         each subject, session, and channel, empty if not analyzed
%     absicoh_within [subjects x sessions x laplace x laplace] - abs(ImCoh) 
%         within hemispheres for each subject and session
%     absicoh_across - same as above but across hemispheres
%     abslagcoh_within [subjects x sessions x laplace x laplace] - abs(LagCoh) 
%         within hemispheres for each subject and session
%     abslagcoh_across - same as above but across hemispheres
%     coh_within [subjects x sessions x laplace x laplace] - coherence
%         within hemispheres for each subject and session
%     coh_across - same as above but across hemispheres
%
% 2. BCI_MI_Laplace_connectivity_results_long.mat - tabular data in the long
%    format for analysis in R (subject ID, session ID, accuracy, 
%    average SNR, measure (ImCoh / LagCoh / Coherence), type (Within / Across),
%    connectivity value)
%%%

analyzed_sessions = zeros(n_subjects, n_sessions);
spec_laplace = zeros(n_subjects, n_sessions, n_laplace, n_freqs);
snr_laplace = NaN(n_subjects, n_sessions, n_laplace);
r2_fooof = NaN(n_subjects, n_sessions, n_laplace);
fooof_fits = cell(n_subjects, n_sessions, n_laplace);
absicoh = NaN(n_subjects, n_sessions, n_laplace, n_laplace);
coh = NaN(n_subjects, n_sessions, n_laplace, n_laplace);
abslagcoh = NaN(n_subjects, n_sessions, n_laplace, n_laplace);

%% Main loop
parfor (subject = 1:n_subjects, num_workers)
% for subject = 1:5
    for session = 1:n_sessions
%     for session = 1
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
        data_bb = double(EEG_bb.data);
        
        % Apply Laplace transform
        data_lap_bb = zeros(n_laplace, size(data_bb, 2), size(data_bb, 3));
        for ch = 1:n_laplace
            ind = laplace{ch, 1};
            ind_neighbors = laplace{ch, 2};
            data_lap_bb(ch, :, :) = data_bb(ind, :, :) - mean(data_bb(ind_neighbors, :, :), 1);            
        end

        % Calculate the spectra
        data_lap_cnt = reshape(data_lap_bb, n_laplace, []);
        spec_lap = pwelch(data_lap_cnt', nfft, noverlap, nfft, EEG.srate)';
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

        % calculate cross-spectra
        tmp_cs = data2cs_bb(data_lap_bb, fres, 1);
        band_freqs = band_freqbins;

        % calculate coherency
        tmp_cohy = squeeze(cs2coh(tmp_cs));

        % calculate derivatives of coherency
        tmp_absicoh = cohy2con(tmp_cohy, 'icoh', 1);
        tmp_coh = cohy2con(tmp_cohy, 'coh', 0);
        tmp_abslagcoh = cohy2con(tmp_cohy, 'lagcoh', 1);

        % average connectivity over frequency band
        [tmp_absicoh_avg, tmp_coh_avg, tmp_abslagcoh_avg] = ...
            aggregate_connectivity({tmp_absicoh, tmp_coh, tmp_abslagcoh}, ...
                                   1, n_laplace, band_freqs);

        absicoh(subject, session, :, :) = tmp_absicoh_avg;
        coh(subject, session, :, :) = tmp_coh_avg;
        abslagcoh(subject, session, :, :) = tmp_abslagcoh_avg;

        analyzed_sessions(subject, session) = 1;
    end
end

disp(sum(analyzed_sessions, [1 2]));

%% Average Connectivity Over Within- and Across-Hemisphere Edges
% Within: C3-CP3 (1-3), C4-CP4 (2-4)
% Across: C3-C4 (1-2), C3-CP4 (1-4), 
%         CP3-C4 (2-3), CP3-CP4 (3-4)
absicoh_within = squeeze(mean(cat(3, squeeze(absicoh(:, :, 1, 3)), ...
    squeeze(absicoh(:, :, 2, 4))), 3));
absicoh_across = squeeze(mean(cat(3, squeeze(absicoh(:, :, 1, 2)), ...
    squeeze(absicoh(:, :, 1, 4)), squeeze(absicoh(:, :, 2, 3)), ...
    squeeze(absicoh(:, :, 3, 4))), 3));

abslagcoh_within = squeeze(mean(cat(3, squeeze(abslagcoh(:, :, 1, 3)), ...
    squeeze(abslagcoh(:, :, 2, 4))), 3));
abslagcoh_across = squeeze(mean(cat(3, squeeze(abslagcoh(:, :, 1, 2)), ...
    squeeze(abslagcoh(:, :, 1, 4)), squeeze(abslagcoh(:, :, 2, 3)), ...
    squeeze(abslagcoh(:, :, 3, 4))), 3));

coh_within = squeeze(mean(cat(3, squeeze(coh(:, :, 1, 3)), ...
    squeeze(coh(:, :, 2, 4))), 3));
coh_across = squeeze(mean(cat(3, squeeze(coh(:, :, 1, 2)), ...
    squeeze(coh(:, :, 1, 4)), squeeze(coh(:, :, 2, 3)), ...
    squeeze(coh(:, :, 3, 4))), 3));

%% Save all results
results_lap_conn = [];
results_lap_conn.analyzed_sessions = analyzed_sessions;

% spectra
results_lap_conn.spec_laplace = spec_laplace;
results_lap_conn.freqs = freqs;

% SNR
results_lap_conn.snr_laplace = snr_laplace;
results_lap_conn.r2_fooof = r2_fooof;
results_lap_conn.fooof_fits = fooof_fits;

% Connectivity
results_lap_conn.absicoh_within = absicoh_within;
results_lap_conn.absicoh_across = absicoh_across;
results_lap_conn.abslagcoh_within = abslagcoh_within;
results_lap_conn.abslagcoh_across = abslagcoh_across;
results_lap_conn.coh_within = coh_within;
results_lap_conn.coh_across = coh_across;

save([savedata 'BCI_MI_Laplace_connectivity.mat'], 'results_lap_conn', '-v7.3');


%% Export to R in the long format
load([savedata 'BCI_MI_task_accuracy.mat']);   % task_accuracy

control_reshape = cell(n_subjects, n_sessions);
accuracy_col = 3;
for subject = 1:n_subjects
    for session = 1:n_sessions
        if results_lap_conn.analyzed_sessions(subject, session) && good_sessions(subject, session)
            control_reshape{subject, session} = ...
                [subject session task_accuracy(subject, session)];
        else
            % mark session with NaN accuracy to be removed later
            control_reshape{subject, session} = ...
                [subject session NaN];
        end
    end
end
metadata = cat(1, control_reshape{:});
snrdata = squeeze(mean(results_lap_conn.snr_laplace, 3));
common_data = repmat(cat(2, metadata, snrdata(:)), 6, 1);
conn_measure = repmat([ ...
    zeros(n_subjects * n_sessions, 1); ...
    ones(n_subjects * n_sessions, 1); ...
    2 * ones(n_subjects * n_sessions, 1)], 2, 1);
conn_type = [zeros(3 * n_subjects * n_sessions, 1); ...
    ones(3 * n_subjects * n_sessions, 1)
];
conn_data = [
    results_lap_conn.absicoh_within(:); 
    results_lap_conn.abslagcoh_within(:); 
    results_lap_conn.coh_within(:); ...
    results_lap_conn.absicoh_across(:); 
    results_lap_conn.abslagcoh_across(:); 
    results_lap_conn.coh_across(:)
];
data = [common_data, conn_measure, conn_type, conn_data];
data(isnan(data(:, accuracy_col)), :) = [];
labels = {'Subject', 'Session', 'Accuracy', 'SNR', ...
    'Measure', 'Type', 'Connectivity'};
measures = {'ImCoh', 'LagCoh', 'Coherence'};
types = {'Within-hemisphere', 'Across-hemisphere'};

save([savedata 'BCI_MI_Laplace_connectivity_results_long.mat'], ...
    'data', 'labels', 'measures', 'types');
