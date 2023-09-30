%%% 
% Fitting CSP to the data during the target presentation interval to
% identify task-relevant voxels
%
% Outputs:
%
% 1. BCI_MI_CSP.mat with the following variables:
%   csp_a_sensor_space - sensor-space CSP patterns
%   csp_a_source_space - source reconstructed CSP patterns
%   csp_cov_matrices - covariance matrices of both classe
%
% 2. BCI_MI_group_CSP.mat with the following variables:
%   csp_w_group, csp_a_group - CSP filters and patterns fitted on averaged
%       covariance matrices
%   Cn - averaged covariance matrices
%
% 3. BCI_MI_CSP_voxel_mask.mat with the following variables:
%   csp_voxel_mask - final mask of task-relevant voxels based on CSP
%   csp_mask - task-relevant voxels for different thresholds
%   csp_roi_mask - task-relevant voxels from ROIs for different thresholds
%   thresholds - considered thresholds
%   num_voxels_roi - number of task-relevant voxels for each ROI
%%%

%% Main loop
A_inverse = {A_eloreta_normal_focal, A_eloreta_normal, A_eloreta_normal_smooth};
n_inverse = numel(A_inverse);
csp_cov_matrices = zeros(n_subjects, n_sessions, n_sensors, n_sensors, n_classes);
csp_a_sensor_space = zeros(n_subjects, n_sessions, n_csp_comps, n_sensors);
csp_a_source_space = zeros(n_subjects, n_sessions, n_inverse, n_csp_comps, n_voxels);

% parfor (subject = 1:n_subjects, num_workers)
for subject = 1:n_subjects
    for session = 1:n_sessions
        desc = ['S' num2str(subject) '/' num2str(session)];
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        
        if isfile([datapath filename])
            % Load the EEG and filter in the mu range
            EEG = pop_loadset('filepath', datapath, 'filename', filename);
            [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, mu_band, downsample, n_mirror_pnts);
            
            % Cut out the target presentation interval, no cursor movement
            EEG_prep = extract_condition_EEG(EEG, EEG_narrow, ...
                {classes_to_analyze 'prep' 'nb'}, events, windows);
            labels = [EEG.etc.trialdata(:).targetnumber];
            epo = eeg2bbci(EEG_prep, labels);
            assert(size(epo.y, 1) == 2, "Expected two classes");
    
            % Save covariance matrices for group-level fit
            csp_cov_matrices(subject, session, :, :, :) = csp_cov(epo);

            % Run CSP for the particular session
            [~, ~, csp_a, ~] = proc_csp(epo, 'SelectFcn', ...
                {@cspselect_equalPerClass, n_csp_comps / 2});
    
            % Convert to source space and normalize
            for iinv = 1:n_inverse
                csp_a_src = csp_a(sa.myinds, :)' * A_inverse{iinv};
                csp_a_src = csp_a_src ./ sum(abs(csp_a_src), 2);
                assert(all(abs(sum(abs(csp_a_src), 2) - 1) < 1e-10), 'source space');
                csp_a_source_space(subject, session, iinv, :, :) = csp_a_src;
            end
    
            % Normalize in sensor space as well
            csp_a = csp_a ./ sum(abs(csp_a), 1);
            assert(all(abs(sum(abs(csp_a), 1) - 1) < 1e-10), 'sensor space');
            csp_a_sensor_space(subject, session, :, :) = csp_a';
    
            % Plot CSP patterns as topomaps in sensor space
            if (csp_make_plots)
                h = figure;
                for c = 1:n_csp_comps
                    subplot(2, n_csp_comps / 2, c); title(['CSP - Comp ', num2str(c)]); 
                    topoplot(csp_a(:, c), EEG.chanlocs); cbar('vert', 0, get(gca, 'clim'));
                end
                sgtitle(desc);
                exportgraphics(h, [savepath_csp 'S' num2str(subject, '%02d') '_Session_' num2str(session, '%02d') '.png']);
                if (run_cmd)
                    close(h);
                end
            end
        end
    end
end

%% Save the results
save([savedata 'BCI_MI_CSP.mat'], ...
    'csp_a_sensor_space', 'csp_a_source_space', 'csp_cov_matrices');

%% Normalize covariance matrices and average
Cn = zeros(n_sensors, n_sensors, n_classes);
for subject = 1:n_subjects
    for session = 1:n_sessions
        if good_sessions(subject, session) > 0
            % normalize trace of the mean covariance matrix before averaging
            % remove differences in total power between subjects and
            % sessions, but keep differences in power between classes
            C = squeeze(csp_cov_matrices(subject, session, :, :, :));
            C_avg = mean(C, 3);
            scale = 1 / trace(C_avg);
%             scale = 1 / norm(C_avg, 'fro');
            Cn = Cn + scale .* C;
            if (trace(C_avg) == 0)
                warning([num2str(subject) '/' num2str(session)]);
            end
        end
    end
end
Cn = Cn / sum(good_sessions, [1 2]);

%% Run CSP analysis for normalized matrices - code was copied from BBCI
[csp_w_group, csp_a_group, ~] = csp_fit_cov(Cn, n_csp_comps / 2);

h = figure;
for c = 1:n_csp_comps
    subplot(2, n_csp_comps / 2, c); title(['CSP - Comp ', num2str(c)]); 
    topoplot(csp_w_group(:, c), all_chanlocs); cbar('vert', 0, get(gca, 'clim'));
end
exportgraphics(h, [savepath_group_csp 'group_csp_filters.png']);

h = figure;
for c = 1:n_csp_comps
    subplot(2, n_csp_comps / 2, c); title(['CSP - Comp ', num2str(c)]); 
    topoplot(csp_a_group(:, c), all_chanlocs); cbar('vert', 0, get(gca, 'clim'));
end
exportgraphics(h, [savepath_group_csp 'group_csp_topomaps.png']);

% Source reconstruction
csp_rec = csp_a_group(sa.myinds, :)' * A_eloreta_normal;
lim = max(abs(csp_rec), [], 2);
h = figure;
ax = subplot(121);
allplots_cortex_subplots(sa, cort5K2full(abs(csp_rec(1, :)), sa), [0 lim(1)], cm17a, 'CSP-R', 1, 'views', -5, 'ax', ax);
ax = subplot(122);
allplots_cortex_subplots(sa, cort5K2full(abs(csp_rec(end, :)), sa), [0 lim(end)], cm17a, 'CSP-L', 1, 'views', -5, 'ax', ax);
exportgraphics(h, [savepath_group_csp 'group_csp_source_rec.png']);

%% Save the results
save([savedata 'BCI_MI_group_CSP.mat'], 'csp_w_group', 'csp_a_group', 'Cn');

%% Plot voxel mask from CSP patterns on the cortex in source space
thresholds = [80 90 95 97.5 99];
n_thresh = numel(thresholds);
h = figure('Position', [10, 10, 917, 437]);
csp_mask = zeros(n_thresh, n_voxels);
csp_roi_mask = zeros(n_thresh, n_voxels);

% Highlight voxels from the target ROIs
roi_mask = zeros(1, n_voxels);
for i_roi = 1:numel(ROI_inds)
    roi_ind = ROI_inds(i_roi);
    
    voxels_roi = get_voxels_roi(sa, roi_ind, []);
    roi_mask(voxels_roi) = 1;
end

for ithresh = 1:n_thresh
    % Select voxels based on the CSP patterns for left and right hand
    % using a threshold
    thresh = thresholds(ithresh);
    prcthresh = prctile(abs(csp_rec), thresh, 2);
    mask = abs(csp_rec) > prcthresh;

    % Use the first and the last CSP patterns to create the mask
    csp_mask(ithresh, mask(1, :)) = 1;
    csp_mask(ithresh, mask(end, :)) = 1;

    % Use only voxels from the sensorimotor ROIs
    csp_roi_mask(ithresh, :) = csp_mask(ithresh, :) & roi_mask;
    
    ax = subplot(2, n_thresh, ithresh);
    allplots_cortex_subplots(sa, cort5K2full(csp_mask(ithresh, :), sa), [0 1], cm17a, 'mask', 1, 'views', 5, 'ax', ax);

    ax = subplot(2, n_thresh, ithresh + n_thresh);
    allplots_cortex_subplots(sa, cort5K2full(csp_roi_mask(ithresh, :), sa), [0 1], cm17a, 'mask', 1, 'views', 5, 'ax', ax);
end

sgtitle({'top: CSP mask, bottom: CSP mask reduced to ROIs', ...
    'left->right: different threshold values - 80 | 90 | 95 | 97.5 | 99'});
exportgraphics(h, [savepath_group_csp 'group_csp_source_thresholded.png']);
if (run_cmd)
    close(h);
end

csp_voxel_mask = squeeze(csp_roi_mask(thresholds == csp_source_threshold, :))';
assert(size(csp_voxel_mask, 1) == n_voxels); % ensure that voxels are the 1st dimension

%% Check which ROIs do selected sources belong to
n_rois_HO = numel(sa.HO_labels) - 1; % ignore Subcortical
num_voxels_roi = zeros(n_rois_HO, 1);
csp_mask_sel = csp_mask(thresholds == csp_source_threshold, :)';
for roi_ind = 1:n_rois_HO
    num_voxels_roi(roi_ind) = numel(get_voxels_roi(sa, roi_ind, csp_mask_sel));
end

%% Save the results
save([savedata 'BCI_MI_CSP_voxel_mask.mat'], 'csp_voxel_mask', ...
    'csp_mask', 'csp_roi_mask', 'thresholds', 'num_voxels_roi');
