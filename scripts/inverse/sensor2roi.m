function [roi_data, ev, w] = sensor2roi(sensor_data, sa, A_inv, agg_method, agg_params)
% SENSOR2ROI Extract time series of activity from ROIs
%
% Parameters:
%   sensor_data [channels x samples (x trials)] - epoched or continuous
%       EEG data
%   sa - source analysis structure
%   A_inv [channels x voxels] - inverse operator
%   agg_method - method for aggregating voxels within the ROI ('svd' or
%       'avg-flip')
%   agg_params - additional parameters for the aggregation of voxels
%
% Returns:
%   roi_data [samples x rois x components] - extracted ROI time series
%   ev [rois x components] - variance explained by the SVD components
%   w {rois x 1} - weights of the SVD components
    [~, n_pnts] = size(sensor_data);

    % Default settings
    n_comps = 1;
    n_rois = numel(sa.HO_labels);
    roi_inds = 1:n_rois;
    voxel_mask = [];
    signflip = [];
    
    if isfield(agg_params, 'n_comps')
        n_comps = agg_params.n_comps;
    end
    if isfield(agg_params, 'roi_inds')
        roi_inds = agg_params.roi_inds;
        n_rois = numel(roi_inds);
    end
    if isfield(agg_params, 'voxel_mask')
        voxel_mask = agg_params.voxel_mask;
    end

    roi_data = zeros(n_pnts, n_rois, n_comps);
    ev = zeros(n_rois, n_comps);
    w = cell(n_rois, 1);
    for i_roi = 1:n_rois
        roi_ind = roi_inds(i_roi);

        voxels_roi = get_voxels_roi(sa, roi_ind, voxel_mask);
        voxel_data = sensor2voxel(sensor_data, sa.myinds, A_inv, voxels_roi);
        
        if strcmp(agg_method, 'avg-flip')
            signflip = agg_params.signflip(voxels_roi);
        end
        
        [roi_data(:, i_roi, :), ev(i_roi, :), w{i_roi}] = ...
            voxel2roi(voxel_data, agg_method, n_comps, signflip);
        
        fprintf('.');
    end
    fprintf('\n');
end