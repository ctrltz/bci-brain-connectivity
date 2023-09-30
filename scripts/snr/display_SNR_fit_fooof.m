function [h] = display_SNR_fit_fooof(fooof_results, snr, fpeak, peak_pos, target_band, ax)
% DISPLAY_SNR_FIT_FOOOF Plot FOOOF fit as a sanity check
%
% Parameters:
%   fooof_results - fitted FOOOF model
%   snr - estimated SNR
%   fpeak - peak frequency or [] if not applicable
%   peak_pos - power of PSD and 1/f at the peak frequency or []
%   target_band - frequency band for the estimation of SNR
%   ax - axes to plot at or [] (new figure will be created if empty)
%
% Returns: 
%   h - figure handle if one was created

    fit_freqs = fooof_results.freqs;
    psd_orig = fooof_results.power_spectrum;
    psd_fit = fooof_results.fooofed_spectrum;
    ap_fit = fooof_results.ap_fit;
    r2 = fooof_results.r_squared;

    h = display_SNR_fit(fit_freqs, psd_orig, psd_fit, ap_fit, ...
        snr, r2, fpeak, peak_pos, target_band, ax);
end