%%% 
% Pre-compute LCMV inverse operators (separate ones for each subject and
% session, using pre-stimulus data from both classes)
%
% Outputs:
% 1. BCI_MI_LCMV.mat with the following variables:
%
%   cov_bb, cov_nb [subjects x sessions x channels x channels] - covariance
%       matrices for broad- and narrowband data
%   A_lcmv_bb, A_lcmv_nb {subjects x sessions}, each cell [channels x
%       voxels] - inverse operators for broad- and narrowband data
%%%

cov_bb = zeros(n_subjects, n_sessions, numel(sa.myinds), numel(sa.myinds));
cov_nb = zeros(n_subjects, n_sessions, numel(sa.myinds), numel(sa.myinds));
A_lcmv_bb = cell(n_subjects, n_sessions);
A_lcmv_nb = cell(n_subjects, n_sessions);

parfor (subject = 1:n_subjects, num_workers)
% for subject = 1:n_subjects
% for subject = [54 43 60]
    for session = 1:n_sessions
%     for session = 1
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];

        if ~isfile([datapath filename])
            continue;
        end

        EEG = pop_loadset('filepath', datapath, 'filename', filename);
        [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
        
        % Calculate covariance matrices using broad- or narrow-band EEG data
        EEG_bb = extract_condition_EEG(EEG, EEG_narrow, {classes_to_analyze, 'rest', 'bb'}, events, windows);
        EEG_nb = extract_condition_EEG(EEG, EEG_narrow, {classes_to_analyze, 'rest', 'nb'}, events, windows);
        data_bb = double(EEG_bb.data);
        data_nb = double(EEG_nb.data);

        C_bb_rest = cov_epo(data_bb, sa.myinds);
        cov_bb(subject, session, :, :) = C_bb_rest;
        C_nb_rest = cov_epo(data_nb, sa.myinds);
        cov_nb(subject, session, :, :) = C_nb_rest;

        % Fit LCMV
        A_lcmv_bb{subject, session} = prepare_LCMV_inverse_operator(C_bb_rest, L_normal, sa.myinds, lambda);
        A_lcmv_nb{subject, session} = prepare_LCMV_inverse_operator(C_nb_rest, L_normal, sa.myinds, lambda);

        fprintf('.');
    end
    fprintf('\n');
end

save("precomputed/BCI_MI_LCMV.mat", "cov_bb", "cov_nb", "A_lcmv_bb", "A_lcmv_nb", '-v7.3');
