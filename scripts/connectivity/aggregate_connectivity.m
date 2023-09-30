function varargout = aggregate_connectivity(data, n_comps, n_rois, band_freqs)
% AGGREGATE_CONNECTIVITY Average connectivity estimates over SVD components
% and frequency bins
%
% Parameters:
%   data - cell array with connectivity matrices
%   n_comps - number of SVD components per ROI
%   n_rois - number of ROIs
%   band_freqs - frequency bins to be aggregated
%
% Returns:
%   aggregated connectivity matrices, the number of outputs is dynamic

    n_outputs = nargout;
    
    varargout = cell(1, n_outputs);
    for idx = 1:n_outputs
        conn = data{idx};
        
        % re-arrange and average across SVD components
        if n_comps > 1
            conn = reshape(conn, [], n_rois, n_comps, n_rois, n_comps);
            conn = squeeze(mean(mean(conn, 3), 5));
        end
        
        % average across frequency band if not narrow-band
        if ~isempty(band_freqs)   % bb
            conn = squeeze(mean(conn(band_freqs, :, :), 1));
        else                      % nb
            conn = squeeze(conn);
        end

        varargout{idx} = conn;
    end
end