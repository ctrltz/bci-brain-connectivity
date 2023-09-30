function [snr, fpeak, dispinfo] = calculate_SNR_fooof(fooof_results, methods, target_band)
% CALCULATE_SNR_FOOOF Calculate SNR using 1/f estimate from FOOOF
%
% Parameters:
%   fooof_results - fitted FOOOF model
%   methods - ways to estimate SNR, can contain one or more of these
%       values:
%       'maxdiff' - maximal difference between PSD and 1/f in log-scale
%       'meandiff' - mean difference between PSD and 1/f in log-scale
%       'ratio' - ratio of PSD and 1/f power (in linear scale)
%   target_band - frequency band for the estimation of SNR
%
% Results:
%   snr [methods x 1] - estimated SNR
%   fpeak [methods x 1] - frequency that corresponds to the maximal SNR
%   snr_pos [method x 2] - power of PSD and 1/f at fpeak (for display
%       purposes)

    n_methods = numel(methods);

    % Unpack FOOOF results    
    fit_freqs = fooof_results.freqs;
    psd_orig = fooof_results.power_spectrum;
    ap_fit = fooof_results.ap_fit;
    
    % Calculate the SNR
    snr = zeros(n_methods, 1);
    fpeak = zeros(n_methods, 1);
    snr_pos = zeros(n_methods, 2);
    band_freqs = fit_freqs >= target_band(1) & fit_freqs <= target_band(2);
    for i = 1:n_methods
        switch (methods{i})
            case 'maxdiff'
                [snr(i), imax] = max(psd_orig(band_freqs) - ap_fit(band_freqs));

                band_freqs_idx = find(band_freqs);
                peak_freq_idx = band_freqs_idx(imax);
                fpeak(i) = fit_freqs(peak_freq_idx);

                snr_pos(i, :) = [ap_fit(peak_freq_idx) psd_orig(peak_freq_idx)];
            case 'meandiff' % averaging log-scaled values might not make a lot of sense
                snr(i) = mean(psd_orig(band_freqs) - ap_fit(band_freqs));
            case 'ratio' % go back to the linear scale and get power ratio
                ptot = mean(10 .^ psd_orig(band_freqs));
                pnoise = mean(10 .^ ap_fit(band_freqs));
                snr(i) = ptot / pnoise;
            otherwise
                error(['unknown method ' method]);
        end
    end

    % Store information for displaying the fit
    dispinfo.fooof_results = fooof_results;
    dispinfo.snr = snr;
    dispinfo.fpeak = fpeak;
    dispinfo.smr_pos = snr_pos;
    dispinfo.target_band = target_band;
end

