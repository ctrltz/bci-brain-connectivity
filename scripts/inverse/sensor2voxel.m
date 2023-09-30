function [voxel_data] = sensor2voxel(sensor_data, smask, A_inv, vmask)
% SENSOR2VOXEL  Obtain source reconstructed time courses of activity
%
% Parameters:
%   sensor_data - sensor-space data, either continuous [channels x samples]
%     or epoched [channels x samples x epochs]
%   smask - indices of channels to use for source reconstruction (should
%     match the provided inverse operator)
%   A_inv - inverse operator for fixed [channels x sources] or free 
%     [channels x sources x 3] dipole orientation 
%   vmask - indices of sources to reconstruct (useful for reducing 
%     amount of computations when analyzing particular ROIs)
%
% Returns:
%   voxel_data - source reconstructed data, shape is determined according to the input
%   arguments

    n_voxels = numel(vmask);
    one_dim_per_voxel = ismatrix(A_inv);
    is_data_continuous = ismatrix(sensor_data);

    if is_data_continuous
        if one_dim_per_voxel            
            voxel_data = tprod(double(sensor_data(smask, :)), [-1 1], A_inv(:, vmask), [-1 2]);
        else
            voxel_data = tprod(double(sensor_data(smask, :)), [-1 1], A_inv(:, vmask, :), [-1 2 3]);
        end
    else
        if one_dim_per_voxel            
            voxel_data = reshape(tprod(double(sensor_data(smask, :, :)), [-1 1 2], A_inv(:, vmask), [-1 3]), [], n_voxels);
        else
            voxel_data = reshape(tprod(double(sensor_data(smask, :, :)), [-1 1 2], A_inv(:, vmask, :), [-1 3 4]), [], n_voxels, 3);
        end
    end
end