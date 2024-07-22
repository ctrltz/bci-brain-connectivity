%%%
% Estimate within- and across-hemisphere phase synchronization for all
% pipelines in the multiverse analysis (Fourier and Hilbert estimation for
% broadband and narrowband data, respectively)
%
% Outputs:
%
% 1. BCI_MI_connectivity_multiverse_BB_NB.mat - results_conn structure with the
% following fields:
%
%     analyzed_sessions [subjects x sessions] - 1 if the session was
%     analyzed, 0 otherwise
%
%     cohy, absicoh, coh, abslagcoh [pipelines x subjects x sessions x rois
%     x rois] - mu band connectivity values for all pairs of ROIs
%
%     absicoh_freq, coh_freq, abslagcoh_freq [pipelines x subjects x
%     sessions x rois x rois x freqs] - frequency-resolved connectivity
%     values for all pairs of ROIs
%
%     absicoh_within, absicoh_across, coh_within, coh_across,
%     abslagcoh_within, abslagcoh_across [pipelines x subjects x sessions]
%     - mu band connectivity values averaged over within- and
%     across-hemisphere connections
%
%%%

analyzed_sessions = zeros(n_subjects, n_sessions);

% Coherency & its derivatives
cohy = complex(zeros(n_multiverse, n_subjects, n_sessions, n_rois, n_rois));
absicoh = zeros(n_multiverse, n_subjects, n_sessions, n_rois, n_rois);
coh = zeros(n_multiverse, n_subjects, n_sessions, n_rois, n_rois);
abslagcoh = zeros(n_multiverse, n_subjects, n_sessions, n_rois, n_rois);

% Frequency-resolved connectivity for plotting the spectra
absicoh_freq = zeros(n_multiverse, n_subjects, n_sessions, n_freqs, n_rois, n_rois);
coh_freq = zeros(n_multiverse, n_subjects, n_sessions, n_freqs, n_rois, n_rois);
abslagcoh_freq = zeros(n_multiverse, n_subjects, n_sessions, n_freqs, n_rois, n_rois);

%% Main loop
parfor (subject = 1:n_subjects, num_workers)
% for subject = 5
    for session = 1:n_sessions
%     for session = 1
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        if ~exist([datapath filename], 'file')
            warning("Skipping subject %d session %d - file not found\n", subject, session);
            continue
        end
        
        %% Prepare EEG
        EEG = pop_loadset('filepath', datapath, 'filename', filename);
        [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
        
        %% Iterate over pipelines
        for i_pipeline = 1:n_multiverse
            pipeline = multiverse{i_pipeline};
            
            disp({subject session multiverse_labels{i_pipeline}});

            %% Extract data for condition
            EEG_condition = extract_condition_EEG(EEG, EEG_narrow, ...
                {classes_to_analyze, 'rest', pipeline.band}, events, windows);
            data_epo = double(EEG_condition.data);    

            %% Prepare the inverse operator
            switch (pipeline.inv_method)
                case 'eLoreta-normal'
                    A_inv = A_eloreta_normal;
                case 'LCMV-normal'
                    if strcmp(pipeline.band, 'bb')
                        A_inv = A_lcmv_bb{subject, session};
                    elseif strcmp(pipeline.band, 'nb')
                        A_inv = A_lcmv_nb{subject, session};
                    else
                        error('Bad value of band')
                    end
                otherwise
                    error('bad inverse method')
            end

            %% Connectivity Calculation
            % Get the data for all combinations of parameters
            [~, ~, n_trials] = size(data_epo);
            agg = pipeline.agg;
            roi_data = sensor2roi(data_epo, sa, A_inv, agg.method, agg);
            roi_data = roi_data(:, :);   % merge roi & comp axes
            % NOTE: in 3SVD case the order is the following (first digit is
            % the index of the ROI, the index of the SVD component)
            % 11 21 31 41 12 22 32 42 13 23 33 43

            % Split data into epochs for connectivity
            n_chans = size(roi_data, 2);
            roi_data_ep = reshape(roi_data', n_chans, [], n_trials);

            % calculate cross-spectra
            if strcmp(pipeline.band, 'bb')
                tmp_cs = data2cs_bb(roi_data_ep, fres, 1);                            
                band_freqs = band_freqbins;
            else
                tmp_cs = data2cs_nb(roi_data_ep, 1);                            
                band_freqs = [];
            end

            % calculate coherency
            tmp_cohy = cs2coh(tmp_cs);

            % calculate derivatives of coherency
            tmp_absicoh = cohy2con(tmp_cohy, 'icoh', 1);
            tmp_coh = cohy2con(tmp_cohy, 'coh', 0);
            tmp_abslagcoh = cohy2con(tmp_cohy, 'lagcoh', 1);

            % x + iy -> abs(x) + i * abs(y) for getting the
            % effective difference to 0/pi
            tmp_cohy = abs(real(tmp_cohy)) + 1i * abs(imag(tmp_cohy));

            % average connectivity over frequency band and SVD
            % components
            [tmp_absicoh_avg, tmp_coh_avg, tmp_abslagcoh_avg, tmp_cohy_avg] = ...
                aggregate_connectivity({tmp_absicoh, tmp_coh, tmp_abslagcoh, tmp_cohy}, ...
                                       agg.n_comps, n_rois, band_freqs);

            cohy(i_pipeline, subject, session, :, :) = squeeze(tmp_cohy_avg);
            absicoh(i_pipeline, subject, session, :, :) = tmp_absicoh_avg;
            coh(i_pipeline, subject, session, :, :) = tmp_coh_avg;
            abslagcoh(i_pipeline, subject, session, :, :) = tmp_abslagcoh_avg;

            % average over SVD components but keep all the
            % frequencies to get connectivity spectra for
            % broadband pipelines
            if strcmp(pipeline.band, 'bb')
                [tmp_absicoh_freq, tmp_coh_freq, tmp_abslagcoh_freq] = ...
                    aggregate_connectivity({tmp_absicoh, tmp_coh, tmp_abslagcoh}, ...
                                           agg.n_comps, n_rois, []);

                absicoh_freq(i_pipeline, subject, session, :, :, :) = tmp_absicoh_freq;
                coh_freq(i_pipeline, subject, session, :, :, :) = tmp_coh_freq;
                abslagcoh_freq(i_pipeline, subject, session, :, :, :) = tmp_abslagcoh_freq;
            end
        end
        
        analyzed_sessions(subject, session) = 1;
    end
end
disp(sum(analyzed_sessions, [1 2]));

%% Average Connectivity Over Within- and Across-Hemisphere Edges
% Within: preL-postL (1-2), preR-postR (3-4)
% Across: preL-preR (1-3), preL-postR (1-4), 
%         postL-preR (2-3), postL-postR (2-4)
absicoh_within = squeeze(mean(cat(4, squeeze(absicoh(:, :, :, 1, 2)), ...
    squeeze(absicoh(:, :, :, 3, 4))), 4));
absicoh_across = squeeze(mean(cat(4, squeeze(absicoh(:, :, :, 1, 3)), ...
    squeeze(absicoh(:, :, :, 1, 4)), squeeze(absicoh(:, :, :, 2, 3)), ...
    squeeze(absicoh(:, :, :, 2, 4))), 4));

abslagcoh_within = squeeze(mean(cat(4, squeeze(abslagcoh(:, :, :, 1, 2)), ...
    squeeze(abslagcoh(:, :, :, 3, 4))), 4));
abslagcoh_across = squeeze(mean(cat(4, squeeze(abslagcoh(:, :, :, 1, 3)), ...
    squeeze(abslagcoh(:, :, :, 1, 4)), squeeze(abslagcoh(:, :, :, 2, 3)), ...
    squeeze(abslagcoh(:, :, :, 2, 4))), 4));

coh_within = squeeze(mean(cat(4, squeeze(coh(:, :, :, 1, 2)), ...
    squeeze(coh(:, :, :, 3, 4))), 4));
coh_across = squeeze(mean(cat(4, squeeze(coh(:, :, :, 1, 3)), ...
    squeeze(coh(:, :, :, 1, 4)), squeeze(coh(:, :, :, 2, 3)), ...
    squeeze(coh(:, :, :, 2, 4))), 4));

%% Save all results
results_conn = [];
results_conn.analyzed_sessions = analyzed_sessions;

results_conn.cohy = cohy;
results_conn.absicoh = absicoh;
results_conn.coh = coh;
results_conn.abslagcoh = abslagcoh;

results_conn.absicoh_freq = absicoh_freq;
results_conn.coh_freq = coh_freq;
results_conn.abslagcoh_freq = abslagcoh_freq;

results_conn.absicoh_within = absicoh_within;
results_conn.absicoh_across = absicoh_across;
results_conn.abslagcoh_within = abslagcoh_within;
results_conn.abslagcoh_across = abslagcoh_across;
results_conn.coh_within = coh_within;
results_conn.coh_across = coh_across;

save([savedata 'BCI_MI_connectivity_multiverse_BB_NB.mat'], 'results_conn', '-v7.3');
