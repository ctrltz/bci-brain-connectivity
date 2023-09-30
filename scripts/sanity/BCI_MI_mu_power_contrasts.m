%%%
% Mu power contrasts between imaginary movements as a sanity check (ERD/ERS
% in sensorimotor areas is expected)
%
% Outputs:
%
% 1. BCI_MI_mu_power_contrasts.mat - results_mu structure with the
% following fields:
%
%   analyzed_sessions [subjects x sessions] - 1 if session was analyzed, 0
%       otherwise
%   spec_sensors [subjects x sessions x periods x classes x sensors x
%       freqs] - channel-wise spectra for each time period and class
%   freqs [freqs x 1] - frequencies corresponding to bins
%   mu_pwr_voxels [subjects x sessions x periods x classes x 2 x voxels] -
%       source space estimates of mu power for each voxel (eLORETA / LCMV)
%
%%%

analyzed_sessions = zeros(n_subjects, n_sessions);

% Spectra
spec_sensors = zeros(n_subjects, n_sessions, n_periods, n_classes, n_sensors, n_freqs);
mu_pwr_voxels = zeros(n_subjects, n_sessions, n_periods, n_classes, n_inverse, n_voxels);

%% Main loop
parfor (subject = 1:n_subjects, num_workers)
% for subject = [5 9 60] 
% for subject = 1:n_subjects
    for session = 1:n_sessions
%     for session = 1
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        if ~exist([datapath filename], 'file')
            warning("Skipping subject %d session %d - file not found\n", subject, session);
            continue
        end
        
        %% Prepare EEG
        EEG = pop_loadset('filepath', datapath, 'filename', filename);
        % narrowband EEG is not actually used for calculating the spectra
        [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
        
        for i_class = 1:n_classes
            for i_period = 1:n_periods
                class = classes_to_analyze(i_class);
                period = periods{i_period};
                    
                disp({subject session class period});

                %% Split data into conditions
                EEG_condition = extract_condition_EEG(EEG, EEG_narrow, {class, period, 'bb'}, events, windows);
                data_epo = double(EEG_condition.data);
                data_cnt = reshape(data_epo, n_sensors, []);

                %% Sensor space
                spec_sensors(subject, session, i_period, i_class, :, :) = pwelch(data_cnt', nfft, noverlap, nfft, EEG.srate)';
                
                %% Source space (voxels)
                for i_inv = 1:n_inverse
                    switch (inv_methods{i_inv})
                        case 'eLoreta-normal'
                            A_inv = A_eloreta_normal;
                        case 'LCMV-normal'
                            A_lcmv = A_lcmv_bb{subject, session};
                        otherwise
                            error('bad inverse method')
                    end

                    source_voxel_data = sensor2voxel(data_epo, sa.myinds, A_inv, 1:n_voxels);

                    spec_voxels = pwelch(source_voxel_data, nfft, noverlap, nfft, EEG.srate);
                    mu_pwr_voxels(subject, session, i_period, i_class, i_inv, :) = squeeze(mean(spec_voxels(mu_band_bins, :), 1));
                end
                   
            end
        end

        analyzed_sessions(subject, session) = 1;
    end
end

disp(sum(analyzed_sessions, [1 2]));

%% Save all results
results_mu = [];
results_mu.analyzed_sessions = analyzed_sessions;
results_mu.spec_sensors = spec_sensors;
results_mu.freqs = freqs;
results_mu.mu_pwr_voxels = mu_pwr_voxels;

save([savedata 'BCI_MI_mu_power_contrasts.mat'], 'results_mu', '-v7.3');