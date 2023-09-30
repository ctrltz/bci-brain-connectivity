function [C] = csp_cov(epo)
% CSP_COV Calculate covariance matrices of two classes for fitting CSP
% (effectively the first part of proc_csp from the BBCI toolbox)
%
% Parameters:
%   epo - epoched data stored in the BBCI format
%
% Returns:
%   C [channels x channels x classes] - covariance matrices

    nChans = size(epo.x, 2);
    nClasses = size(epo.y, 1);
    C = zeros(nChans, nChans, nClasses);
    for k = 1:nClasses
        idx = find(epo.y(k, :));
        X = permute(epo.x(:, :, idx), [1 3 2]);
        X = reshape(X, [], nChans);
        C(:, :, k) = cov(X);
    end
end