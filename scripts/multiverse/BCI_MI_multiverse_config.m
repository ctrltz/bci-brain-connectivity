%%% 
% Configuration and the list of considered options for the multiverse 
% analysis
%%% 

inv_methods = {'eLoreta-normal', 'LCMV-normal'};
n_inverse = numel(inv_methods);

% Methods for extraction of ROI time series 
roi_methods = {'1SVD', '3SVD', 'AVG-flip'};
n_roi_methods = numel(roi_methods);
n_roi_defs = 2;  % mask, no mask
n_roi_comps = 3; % max 3 comps per ROI
roi_agg = {
    struct('label', '1SVD', 'method', 'svd', 'n_comps', 1, 'roi_inds', ROI_inds), ...
    struct('label', '3SVD', 'method', 'svd', 'n_comps', 3, 'roi_inds', ROI_inds), ...
    struct('label', 'AVG-flip', 'method', 'avg-flip', 'n_comps', 1, 'roi_inds', ROI_inds, 'signflip', signflip), ...
    struct('label', 'mask/1SVD', 'method', 'svd', 'n_comps', 1, 'roi_inds', ROI_inds, 'voxel_mask', csp_voxel_mask), ...
    struct('label', 'mask/3SVD', 'method', 'svd', 'n_comps', 3, 'roi_inds', ROI_inds, 'voxel_mask', csp_voxel_mask), ...
    struct('label', 'mask/AVG-flip', 'method', 'avg-flip', 'n_comps', 1, 'roi_inds', ROI_inds, 'voxel_mask', csp_voxel_mask, 'signflip', signflip)
};
n_roi_agg = numel(roi_agg);
roi_agg_labels = cellfun(@(x) x.label, roi_agg, 'UniformOutput', false);
assert(n_roi_agg == n_roi_methods * n_roi_defs);

% Options for the multiverse analysis
n_multiverse = n_inverse * n_bands * n_roi_agg;
multiverse = cell(n_multiverse, 1);
multiverse_labels = cell(n_multiverse, 1);
multiverse_data = zeros(n_multiverse, 5);
k = 1;
fprintf('Multiverse:\n');
for i_inv = 1:n_inverse
    for i_band = 1:n_bands
        for i_agg = 1:n_roi_agg
            multiverse{k} = struct('inv_method', inv_methods{i_inv}, 'band', bands{i_band}, 'agg', roi_agg{i_agg});
            multiverse_labels{k} = [multiverse{k}.inv_method '/' multiverse{k}.band '/' multiverse{k}.agg.label];
            agg = roi_agg{i_agg};
            % inverse method, band, aggregation method, no mask (0) / mask
            % (1), extraction method: 1SVD (1) / 3SVD (2) / AVG-flip (3)
            multiverse_data(k, :) = [i_inv i_band i_agg (i_agg > 3) 1+mod(i_agg-1, 3)];
            fprintf('%d: %s\n', k, multiverse_labels{k});
            k = k + 1;
        end
    end
end