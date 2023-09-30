function [cs] = data2cs_bb(data, fres, segave)
%DATA2CS_BB Compute cross-spectra from epoched data via FFT
%
% Parameters:
%   data - epoched data [channels x samples x epochs]
%   fres - frequency resolution for FFT (number of bins in (0; nyquist])
%   segave - return average CS if 1, epoch-wise CS if not
%
% Returns:
%   estimated cross-spectra, optionally averaged over epochs
%   segave = 1 -> cs [frequencies x channels x channels]
%   segave = 0 -> cs [frequencies x channels x channels x epochs]
%
    if nargin < 3
        segave = 1;
    end

    [nchan, ndat, nep] = size(data);
    
    maxfreq = fres + 1;
    window = ndat;
    overlap = round(ndat / 2);

    cs = zeros(nchan, nchan, maxfreq, nep);
    
    % as in MVGC
    data = demean(data);
    data = permute(data, [2 1 3]); % transpose row, col (works for single-trial data too)

    for r = 1:nep % works for single-trial too
        cs(:, :, :, r) = cpsd_welch(data(:, :, r), nchan, fres+1, window, overlap);
    end

    cs = cs(:, :, 1:maxfreq, :);
    cs = permute(cs, [3, 1, 2, 4]);
    
    if segave
        cs = squeeze(mean(cs, 4));
    end
end


% copied from MVGC toolbox
function S = cpsd_welch(X,n,h,window,noverlap)

    nfft = 2*(h-1);

    S = complex(zeros(n,n,h));

    for i = 1:n
        S(i,i,:) = pwelch(X(:,i),window,noverlap,nfft);          % auto-spectra
        for j = i+1:n % so we don't compute cross-spectra twice
            S(i,j,:) = cpsd(X(:,i),X(:,j),window,noverlap,nfft); % cross-spectra
            S(j,i,:) = conj(S(i,j,:));
        end
    end
end

