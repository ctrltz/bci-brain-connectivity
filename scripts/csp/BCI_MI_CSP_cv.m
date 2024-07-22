%%%
% Train a CSP + rLDA classifier to distinguish between left- and right-hand
% imaginary movements using data from each time window
%
% Outputs:
% 1. BCI_MI_CSP_cv.mat - accuracy and AUC (6-fold CV) + classifier settings
% 2. BCI_MI_CSP_accuracy.mat - same data in long format for R
%%%

proc.train = {
    {'CSPW', @proc_csp, 'SelectFcn', ...
        {@cspselect_equalPerClass, n_csp_comps / 2}}, ...
    @proc_variance, ...
    @proc_logarithm
};
proc.apply = {
    {@proc_linearDerivation, '$CSPW'}, ...
    @proc_variance, ...
    @proc_logarithm
};

accuracy_cv = {
    NaN(n_subjects, n_sessions), ...
    NaN(n_subjects, n_sessions), ...
    NaN(n_subjects, n_sessions)
};
auc_cv = {
    NaN(n_subjects, n_sessions), ...
    NaN(n_subjects, n_sessions), ...
    NaN(n_subjects, n_sessions)
};

for subject = 1:n_subjects
    for session = 1:n_sessions
        desc = ['S' num2str(subject) '/' num2str(session)];
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];

        if isfile([datapath filename])
            % Load the EEG and filter in the mu range
            EEG = pop_loadset('filepath', datapath, 'filename', filename);
            [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
            
            for i_period = 1:n_periods
                disp({subject session periods{i_period}});

                % Cut out the time interval of interest
                EEG_cv = extract_condition_EEG(EEG, EEG_narrow, ...
                    {classes_to_analyze periods{i_period} 'nb'}, events, windows);
                labels = [EEG.etc.trialdata(:).targetnumber];
                epo = eeg2bbci(EEG_cv, labels);
                assert(size(epo.y, 1) == 2, "Expected two classes");
    
                loss_mean = ...
                    crossvalidation(epo, {@train_RLDAshrink, 'Gamma', 0.05}, ...
                                    'SampleFcn', {@sample_KFold, 6}, ...
                                    'LossFcn', {@loss_0_1, @loss_rocArea}, ...
                                    'Proc', proc);
    
                accuracy_cv{i_period}(subject, session) = 1 - loss_mean(1);
                auc_cv{i_period}(subject, session) = 1 - loss_mean(2);
            end
        end
    end
end

%% Save the results
save([savedata 'BCI_MI_CSP_cv.mat'], 'accuracy_cv', 'auc_cv', 'proc');


%% Export to R
% Load the online accuracy values
load([savedata 'BCI_MI_task_accuracy.mat']);   % task_accuracy

control_reshape_subject = repmat((1:n_subjects)', 1, n_sessions);
control_reshape_session = repmat(1:n_sessions, n_subjects, 1);

csp_accuracy_data = [];
for i_period = 1:n_periods
    period_data = [
        control_reshape_subject(:) control_reshape_session(:) ...
        i_period * ones(n_subjects * n_sessions, 1) ...
        accuracy_cv{i_period}(:) auc_cv{i_period}(:) ...
        task_accuracy(:) good_sessions(:)
    ];
    csp_accuracy_data = [csp_accuracy_data; period_data];
end

% Remove missing and excluded sessions
csp_accuracy_data = csp_accuracy_data(csp_accuracy_data(:, 7) > 0, 1:6);
labels = {'Subject', 'Session', 'Period', 'Accuracy', 'AUC', 'Online'};

save([savedata 'BCI_MI_CSP_accuracy.mat'], ...
    'csp_accuracy_data', 'labels', 'periods');