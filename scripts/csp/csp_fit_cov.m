function [W, A, score] = csp_fit_cov(C, n_comps)
% CSP_FIT_COV Fit CSP using provided covariance matrices of the classes
% (effectively the second part of proc_csp from the BBCI toolbox)
%
% Parameters:
%   C [channels x channels x classes] - covariance matrices
%   n_comps - number of CSP components per class to return
%
% Returns:
%   W [channels x comps] - CSP filters
%   A [channels x comps] - CSP patterns
%   score - eigenvalues

    % Whitening
    M = procutil_whiteningMatrix([], 'C', mean(C, 3));
    if (size(M, 2) < size(C, 1))
        warning('Due to dimensionality reduction a maximum of only %d CSP components can be computed', size(M,2))
    end
    
    % Eigenvalue decomposition
    [W, D] = eig(M' * (C(:,:,1) - C(:,:,2)) * M);
    W = M * W; % project filters from whitened space back into original channel space
    [ev, sort_idx] = sort(diag(D), 'ascend');
    D = diag(ev);
    W = W(:,sort_idx);
    
    % Select equal number of components for each class
    idx = cspselect_equalPerClass(ev, W, n_comps);
    W = W(:, idx);
    W = W ./ vecnorm(W, 2, 1);
    score = ev(idx);
    
    C_avg = mean(C, 3);
    A = C_avg * W / (W'*C_avg*W);
end

