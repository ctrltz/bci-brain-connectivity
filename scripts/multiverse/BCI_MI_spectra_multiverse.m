%%%
% Calculate ROI spectra for all pipelines in the multiverse analysis
%
% Output:
%
% 1. BCI_MI_spectra_rois_multiverse.mat - results_spec structure with the 
%    following fields:
%    
%    analyzed_sessions [subjects x sessions] - 1 if session was analyzed, 0
%        otherwise
%    spec_rois {pipelines x subjects x sessions} - ROI spectra for each
%        pipeline, subject, and session, [] if not analyzed
%    ev_rois [pipelines x subjects x sessions x rois x 3] - variance
%        explained by each of the SVD components
%    w_rois {pipelines x subjects x sessions x rois} - weights for each of
%        the SVD components
%    freqs [freqs x 1] - frequencies corresponding to the bins
%%%

tic;

%% Main loop
analyzed_sessions = zeros(n_subjects, n_sessions);

spec_rois = cell(n_multiverse, n_subjects, n_sessions);
ev_rois = zeros(n_multiverse, n_subjects, n_sessions, n_rois, n_roi_comps);
w_rois = cell(n_multiverse, n_subjects, n_sessions, n_rois);

parfor (subject = 1:n_subjects, num_workers)
%     for subject = 1:12 
    for session = 1:n_sessions
%         for session = 1
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        if ~exist([datapath filename], 'file')
            warning("Skipping subject %d session %d - file not found\n", subject, session);
            continue
        end
        
        %% Prepare EEG
        EEG = pop_loadset('filepath', datapath, 'filename', filename);
        % narrowband EEG is not actually used for calculating the spectra
        [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
        
        %% Iterate over pipelines
        for i_pipeline = 1:n_multiverse
            pipeline = multiverse{i_pipeline};
            
            disp({subject session multiverse_labels{i_pipeline}});

            %% Extract pre-stimulus data
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
                    error('Bad inverse method')
            end
            
            %% Power spectral density
            [~, ~, n_trials] = size(data_epo);
            agg = pipeline.agg;
            [roi_data, ev, w] = sensor2roi(data_epo, sa, A_inv, agg.method, agg);
            roi_data = roi_data(:, :);   % merge roi & comp axes
                
            spec_rois{i_pipeline, subject, session} = pwelch(...
                roi_data, nfft, noverlap, nfft, srate)';
            if strcmp(agg.method, 'svd') && agg.n_comps == 3  % save info from 3SVD
                ev_rois(i_pipeline, subject, session, :, :) = ev;
                w_rois(i_pipeline, subject, session, :) = w;
            end
        end
                   
        analyzed_sessions(subject, session) = 1;
    end
end

disp(sum(analyzed_sessions, [1 2]));

%% Save all results
results_spec = [];
results_spec.analyzed_sessions = analyzed_sessions;
results_spec.spec_rois = spec_rois;
results_spec.ev_rois = ev_rois;
results_spec.w_rois = w_rois;
results_spec.freqs = freqs;

save([savedata 'BCI_MI_spectra_rois_multiverse.mat'], 'results_spec', '-v7.3');