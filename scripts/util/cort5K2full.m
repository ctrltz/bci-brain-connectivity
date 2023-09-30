function [data_remapped] = cort5K2full(data, sa)
% CORT5K2FULL Map data from 4502 voxels in cortex5K to cortex75K for 
% visualization
%
% Parameters:
%   data - data to be mapped (4502 values are expected)
%   sa - source analysis structure
%
% Returns:
%   data_remapped - same data, but copied to fill the 75K space

    data = squeeze(data);
    data_remapped = zeros(numel(sa.cortex5K.in_from_cortex75K), 1); 
    assert(numel(data) == numel(sa.voxels_5K_cort));

    % Map 4502 values (in cortex) into 5004 values (5K model)
    data_remapped(sa.cortex5K.in_cort) = data;

    % Map 5004 values (5K model) into 74382 values (75K model)
    data_remapped = data_remapped(sa.cortex5K.in_to_cortex75K_geod);
end

