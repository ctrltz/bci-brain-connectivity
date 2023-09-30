function [h] = display_SNR_fit(freqs, psd, yhat, noise, snr, r2, fpeak, snr_pos, target_band, ax)
% DISPLAY_SNR_FIT Plot 1/f fit and SNR estimate as a sanity check
%
% Parameters:
%   freqs [freqs x 1] - vector of frequencies used for estimation of the spectra
%   psd [freqs x 1] - estimated PSD
%   yhat [freqs x 1] - fitted PSD
%   noise [freqs x 1] - fitted aperiodic component
%   snr - estimated SNR
%   r2 - r2 of the 1/f fit
%   fpeak - peak frequency or [] if not applicable
%   snr_pos - power of PSD and 1/f at the peak frequency or []
%   target_band - frequency band for the estimation of SNR
%   ax - axes to plot at or [] (new figure will be created if empty)
%
% Returns: 
%   h - figure handle if one was created

    % Create figure if needed
    h = [];
    if isempty(ax)
        h = figure; 
        ax = gca;
    end    
    hold on;
    
    % plot spectra and noise fit
    p1 = plot(ax, freqs, psd, 'k', 'LineWidth', 1.5);
    p3 = plot(ax, freqs, noise, 'b', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('log10(PSD)');
    
    % display target band
    if ~isempty(target_band)
        xline(ax, target_band, '--');
    end
    
    % add SMR position display
    if ~isempty(fpeak) && ~isempty(snr_pos)
        line(ax, [fpeak fpeak], snr_pos, 'Color', 'm', 'LineWidth', 1.5);
        text(ax, fpeak, snr_pos(2) + snr * 0.3, {[num2str(fpeak) ' Hz'], [num2str(snr) ' dB']});
    end
    
    % plot PSD fit if available
    if ~isempty(yhat)
        p2 = plot(ax, freqs, yhat, 'r', 'LineWidth', 1.5);
        legend(ax, [p1 p2 p3], {'PSD', 'PSD fit', 'noise fit'}, 'Location', 'northeast');
    else
        legend(ax, [p1 p3], {'PSD', 'noise fit'}, 'Location', 'northeast');
    end
    
    % display the SNR in title
    title(ax, ['SNR = ', num2str(snr, '%.2f'), ' | r^2 = ', num2str(r2, '%.2f')]);
end

