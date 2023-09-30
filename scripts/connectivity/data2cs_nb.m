function [cs] = data2cs_nb(data, segave)
% DATA2CS_NB Compute cross-spectra from epoched data via analytic signal
%
% Parameters:
%   data - epoched data [channels x samples x epochs]
%   segave - return average CS if 1, epoch-wise CS if not
%
% Returns:
%   estimated cross-spectra, optionally averaged over epochs
%   frequency axis is kept for consistency with broadband version
%   segave = 1 -> cs [1 x channels x channels]
%   segave = 0 -> cs [1 x channels x channels x epochs]

    if nargin < 2
        segave = 1;
    end

    [nchan, ~, nep] = size(data);
    cs = zeros(nchan, nchan, nep);

    % as in MVGC
    data = demean(data);
    data = permute(data, [2 1 3]); % transpose row, col (works for single-trial data too)
    
    for r = 1:nep
        dataep = hilbert(squeeze(data(:, :, r)));
        cs(:, :, r) = conj(dataep' * dataep);
    end
    
    if segave
        cs = squeeze(mean(cs, 3));
    end

    % add first dimension (freq) for consistency with broadband case
    cs = reshape(cs, [1 size(cs)]);
end