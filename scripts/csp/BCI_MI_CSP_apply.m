%%% 
% Calculate the spectra of CSP components to show class difference
%
% Outputs:
% 1. BCI_MI_CSP_spectra.mat with the following variables:
%   spec_csp [subjects x sessions x classes x comps x freqs] - spectra of
%       all CSP components for each subject, session, and class
%   spec_csp_avg [classes x comps x freqs] - averaged over all subjects and
%       sessions
%%%

%% Main loop
spec_csp = zeros(n_subjects, n_sessions, n_classes, n_csp_comps, n_freqs);

% parfor (subject = 1:n_subjects, num_workers)
for subject = 1:n_subjects
    for session = 1:n_sessions
        desc = ['S' num2str(subject) '/' num2str(session)];
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        
        if isfile([datapath filename])
            % Load the EEG data
            EEG = pop_loadset('filepath', datapath, 'filename', filename);
            [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
            
            % Cut out the target presentation interval, no cursor movement
            % Use broadband data for calculating the spectra
            EEG_prep = extract_condition_EEG(EEG, EEG_narrow, ...
                {classes_to_analyze 'prep' 'bb'}, events, windows);
            labels = [EEG.etc.trialdata(:).targetnumber];
            
            % Get spectra for both classes
            for c = classes_to_analyze
                disp(c);
                class_trials = find(labels == c);
                EEG_class = pop_select(EEG_prep, 'trial', class_trials);
                X = reshape(double(EEG_class.data), n_sensors, []);
                X_csp = csp_w_group' * X;
                spec_csp(subject, session, c, :, :) = ...
                    pwelch(X_csp', nfft, noverlap, nfft, srate)';
            end
        end
    end
end

%% Calculate average spectra for both classes
good_sessions_flat = reshape(good_sessions, [], 1);
spec_csp_flat = reshape(spec_csp, [], n_classes, n_csp_comps, n_freqs);
spec_csp_avg = squeeze(mean(spec_csp_flat(good_sessions_flat > 0, :, :, :), 1));

%% Save the results
save([savedata 'BCI_MI_CSP_spectra.mat'], 'spec_csp', 'spec_csp_avg');