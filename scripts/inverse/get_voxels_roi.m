function [voxels_roi] = get_voxels_roi(sa, roi_ind, voxel_mask)
% GET_VOXELS_ROI Extract indices of voxels in a given ROI
%
% Parameters:
%   sa - source analysis structure
%   roi_ind - index of the ROI
%   voxel_mask - mask to consider only a subset of voxels
%
% Returns:
%   voxels_roi - indices of voxels that belong to the ROI

    roi_mask = sa.cortex75K.in_HO(sa.voxels_5K_cort) == roi_ind;
    if ~isempty(voxel_mask) 
        roi_mask = roi_mask & voxel_mask;
    end
    voxels_roi = find(roi_mask > 0);
end