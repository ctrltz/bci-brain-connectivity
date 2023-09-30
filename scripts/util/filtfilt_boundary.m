function [filtered_data] = filtfilt_boundary(b, a, data, event, n_mirror_pnts)
% FILTFILT_BOUNDARY Filter the data taking into account boundary events
% (mirroring is performed to reduce edge effects)
%
% Parameters:
%   b, a - parameters of the filter (as for usual filtfilt)
%   data [samples x channels] - data to be filtered (as for usual filtfilt)
%   event - list of events to consider (only 'boundary' ones are used)
%   n_mirror_pnts - number of samples to use for mirroring
%
% Returns:
%   filtered_data [samples x channels] - data after filtering, shape is the
%   same as data

    is_boundary = strcmp({event(:).type}, 'boundary');
    boundary_lat = [event(is_boundary).latency];
    cont_segments = [0 boundary_lat; boundary_lat size(data, 1)];
    n_segments = size(cont_segments, 2);

    fprintf('Filtering %d continuous segments...', n_segments);
    filtered_data = zeros(size(data));
    for i = 1:n_segments
        s_idx = round(cont_segments(1, i)) + 1;
        e_idx = round(cont_segments(2, i));
        
        % Mirror the signal to suppress edge effects due to filtering
        n_mirror_actual = min(n_mirror_pnts, e_idx - s_idx);

        data_prefix = 2 * data(s_idx, :) - flipud(data(s_idx+1:s_idx+n_mirror_actual, :));
        data_suffix = 2 * data(e_idx, :) - flipud(data(e_idx-n_mirror_actual:e_idx-1, :));
        data_mirrored = [data_prefix; data(s_idx:e_idx, :); data_suffix];

        % Filter and use the original part
        data_filtered = filtfilt(b, a, data_mirrored);
        filtered_data(s_idx:e_idx, :) = data_filtered(n_mirror_actual+1:end-n_mirror_actual, :);
    end
    fprintf('Done\n');
end