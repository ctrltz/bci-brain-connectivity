function [A_lcmv] = prepare_LCMV_inverse_operator(C, L, smask, lcmv_reg)
% PREPARE_LCMV_INVERSE_OPERATOR Fit LCMV on the provided covariance matrix
% of the data
%
% Parameters:
%   C [channels x channels] - covariance matrix
%   L [channels x voxels] - lead field matrix (fixed orientations)
%   smask - mask to use only a subset of channels
%   lcmv_reg - regularization parameter for LCMV (0.05 if not provided)
%
% Returns:
%   A_lcmv [channels x voxels] - LCMV inverse operator

    % Regularization constant for LCMV
    if nargin < 4
        lcmv_reg = 0.05;
    end

    n_chans = numel(smask);       
    alpha = lcmv_reg * trace(C) / length(C);
    Cr = C + alpha * eye(n_chans);
    [~, A_lcmv] = lcmv(Cr, L, struct('alpha', 0, 'onedim', 1));
end