function [channels_excluded, trials_excluded, hs] = extract_spectra_outliers_auto(EEG, ...
    f_range, crit, thresh_chan, thresh_ep, sensitive_channels, ...
    exclude_channels, exclude_trials, exclude_sensitive, plot)
% EXTRACT_SPECTRA_OUTLIERS_AUTO Identify bad trials and channels as outliers
% (abs(z) > criterion) based on power in a certain frequency range
% (recursive procedure, decision is made automatically)
%
% Parameters:
%   EEG - data in EEGLAB
%   f_range - frequency ranges for power calculation
%   crit - criterion (currently works only for the first band)
%   thresh_chan - maximal percentage of outlier trials for good channels
%   thresh_ep - maximal percentage of outlier channels for good trials
%   sensitive_channels - epochs are bad if these channels are outliers
%   exclude_channels - exclude channels if 1
%   exclude_trials - exclude trials if 1
%   exclude_sensitive - exclude trials based on sensitive channels if 1
%   plot - plot z-scores on each iteraction if 1
%
% Returns:
%   channels_excluded - bad channels
%   trials_excluded - bad trials
%   hs - handles to all the generated figures
    
    nfft = 2 ^ (nextpow2(EEG.srate) + 1);
    noverlap = nfft / 2;
    freqs = psdfreqvec('npts', nfft, 'Fs', EEG.srate);

    [nf, n_bands] = size(f_range);
    assert(nf == 2);
    
    % store all the figures
    hs = [];
    
    % convert frequencies to bins
    fb_range = zeros(size(f_range));
    for b = 1:n_bands
        fb_range(:, b) = dsearchn(freqs, f_range(:, b));
    end
    
    sensitive_chans = ismember({EEG.chanlocs(:).labels}, sensitive_channels);
    trial_start = EEG.event(ismember({EEG.event(:).type}, {'0'}));
    trial_end = EEG.event(ismember({EEG.event(:).type}, {'3'}));
    assert(numel(trial_start) == numel(trial_end), "Numbers of start and end triggers do not match");
    n_trials = numel(trial_start);
    if n_trials ~= 450
        warning(['Number of trials is not 450, but ', num2str(n_trials)]);
    end
    pwr = zeros(EEG.nbchan, n_trials, n_bands);
    
    for i = 1:n_trials
        slat = floor(trial_start(i).latency);
        % downsampling led to out-of-bound events, 
        elat = ceil(trial_end(i).latency) - 1;
        dataep = EEG.data(:, slat:elat);
        [spec_ep, ~] = pwelch(dataep', nfft, noverlap, nfft, EEG.srate);
        spec_ep = 10*log10(spec_ep);
        
        for b = 1:n_bands
            pwr(:, i, b) = squeeze(mean(spec_ep(fb_range(1, b):fb_range(2, b), :), 1));
        end
    end
    
    channels_excluded = [];
    trials_excluded = [];
    iter = 1;
    stop = false;
    
    while (~stop)
        fprintf('Iteration %d...', iter);
        
        % Zero out channels and epochs that already were excluded
        pwr_iter = pwr;
        pwr_iter(channels_excluded, :, :) = 0;
        pwr_iter(:, trials_excluded, :) = 0;

        zpwr = zscore(pwr_iter, 0, [1 2]);
        zpwr(channels_excluded, :, :) = NaN;
        zpwr(:, trials_excluded, :) = NaN;

        if (plot)
            hc = figure;
            hc.WindowState = 'maximized';
            pause(0.1);
            for b = 1:n_bands
                subplot(n_bands, 1, b);
                imagesc(squeeze(zpwr(:, :, b))); axis xy; colorbar;
                xlabel('Epoch'); ylabel('Channel');
                title([num2str(f_range(1, b)) ' - ' num2str(f_range(2, b)) ' Hz']);
            end
            sgtitle(['Iteration ' num2str(iter)]);
            hs = [hs hc];
        end
        
        % Check if epochs can be excluded
        trial_outliers = squeeze(mean(abs(zpwr(:, :, 1)) > crit, 1, 'omitnan')) > thresh_ep;
        if exclude_trials && any(trial_outliers)
            trial_outliers = find(trial_outliers);
            fprintf('trials to be excluded:\n');
            disp(trial_outliers);
            trials_excluded = [trials_excluded; trial_outliers'];
            iter = iter + 1;
            continue;
        end       
        
        % Check if channels can be excluded
        channel_outliers = squeeze(mean(abs(zpwr(:, :, 1)) > crit, 2, 'omitnan')) > thresh_chan;
        if exclude_channels && any(channel_outliers)
            channel_outliers = find(channel_outliers);
            fprintf('channels to be excluded:\n');
            disp(channel_outliers');
            channels_excluded = [channels_excluded; channel_outliers];
            iter = iter + 1;
            continue;
        end

        % Check if sensitive channels contain artifacts
        if (exclude_sensitive)
            sensitive_outliers = any(squeeze(zpwr(sensitive_chans, :, 1)) > crit, 1);
            if any(sensitive_outliers)
                sensitive_outliers = find(sensitive_outliers);
                fprintf('trials to be excluded due to artifacts in sensitive channels:\n');
                disp(sensitive_outliers);
                trials_excluded = [trials_excluded; sensitive_outliers'];
                iter = iter + 1;
                continue;
            end
        end
        
        % Stop if found nothing to exclude
        fprintf('Nothing needs to be excluded. Stopping\n');
        stop = true;
    end
end

