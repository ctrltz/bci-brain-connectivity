function [roi_data, ev, w] = voxel2roi(voxel_data, agg_method, n_comps, signflip)
% VOXEL2ROI Obtain time courses of source reconstructed activity in a ROI
%
% Parameters:
%   voxel_data - source-space data, either continuous [channels x samples]
%     or epoched [channels x samples x epochs]
%   agg_method ['svd' | 'avg' | 'avg_flip'] - method for aggregation of
%     time series of activity for sources within the ROI
%   n_comps - how many components per ROI to extract (only used when
%     agg_method is 'svd', otherwise will effectively be 1)
%   signflip - vector of signs (-1 or 1) that should be applied to sources
%     within the ROI before averaging (only used when agg_method is
%     'avg_flip')
%
% Returns:
%   roi_data - source reconstructed data, third dimension is flattened if
%     epochs were provided [n_comps x (samples*epochs)]
%   ev - variance explained by each of the kept SVD components (makes sense
%     only if agg_method is 'svd')
%   w - source weights corresponding to the kept SVD components (make sense
%     only if agg_method is 'svd')

    % explained variance and weights are meaningful only for SVD
    ev = 0;
    w = [];

    switch (agg_method)
        case 'svd'
            voxel_data = voxel_data(:, :) - mean(voxel_data(:, :), 1);
            [voxel_svd, sigma, coeff] = svd(voxel_data(:, :), 'econ');

            roi_data = voxel_svd(:, 1:n_comps);
            ev_tmp = diag(sigma .^ 2) ./ sum(diag(sigma .^ 2));
            assert(abs(1 - sum(ev_tmp)) < 1e-10);
            ev = ev_tmp(1:n_comps);
            w = coeff(:, 1:n_comps);
        case 'avg'
            roi_data = mean(voxel_data(:, :), 2);
        case 'avg-flip'
            [n_samples, ~, n_dims] = size(voxel_data);
            flipmask = repmat(signflip', n_samples, 1, n_dims);
            roi_data = mean(voxel_data(:, :) .* flipmask(:, :), 2);
        otherwise
            error(['Method ' agg_method ' is not supported for ROI time series aggregation']);
    end
end