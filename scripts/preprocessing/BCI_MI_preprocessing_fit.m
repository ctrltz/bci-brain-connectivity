function [] = BCI_MI_preprocessing_fit(subject, session, par)
% BCI_MI_PREPROCESSING_FIT Preprocess the data
%
% Parameters:
%   subject, session
%   par - containter for multiple parameters unpacked below

    %% Unpack par
    datapath = par.data.raw;    % path to the raw data
    auxpath = par.preproc.aux;  % path to the preprocessing info
    savedata = par.preproc.ica; % path to the intermediate data
    savepath = par.savepath;    % path to the figures
    ss = par.ss;                % statuses enum
    stage = par.stage;          % statuses for every stage

    %% Prepare file and folder names
    subjprefix_plot = ['S' num2str(subject) '/', num2str(session)];
    subjprefix = ['S' num2str(subject) '_Session_', num2str(session)];
    subjfolder = ['fit/' subjprefix '/'];
    good_chans_prefix = 'good_channels/';
    bad_chans_prefix = 'bad_channels/';
    good_comps_prefix = 'good_components/';
    bad_comps_prefix = 'bad_components/';
    mkdir_if_not_exists([savepath subjfolder]);
    
    %% Load data
    filename = [subjprefix '.mat'];
    if ~exist([datapath filename], 'file')
        warning([filename ' - file not found, skipping']);
        return
    end
    fprintf('Loading data from %s... ', filename);
    ses_data = load([datapath filename]);
    fprintf('Done\n');

    %% Load already available preprocessing info
    preproc_info_filename = [subjprefix '_preproc_info.mat'];
    if exist([auxpath preproc_info_filename], 'file')
        tmp = load([auxpath preproc_info_filename]);
        preproc_info = tmp.preproc_info;
    else
        preproc_info = [];
    end

    %% Stage 1: rejecting bad channels and epochs
    if (stage.bad_channels_epochs == ss.apply)
        fprintf('[%s - stage 1 - apply] rejecting previously identified bad channels and epochs\n', subjprefix_plot);
        EEG = bci2eeglab(ses_data.BCI, 1, par.electrode_file, 1, 1);

        % Down-sample the data for speed
        EEG = pop_resample(EEG, 250); % in Hz
        EEG = eeg_checkset(EEG);

        % Remove previously identified bad trials
        if ~isempty(preproc_info.bad_trials)
            epochs_excluded = trials2epochs(EEG, preproc_info.bad_trials);
            EEG = pop_select(EEG, 'nopoint', epochs_excluded);
            EEG = eeg_checkset(EEG); % check data for inconsistencies
        end

        % Remove previously identified bad channels
        if ~isempty(preproc_info.bad_channels)
            EEG = pop_select(EEG, 'nochannel', preproc_info.bad_channels);
            EEG = eeg_checkset(EEG); % check data for inconsistencies
        end

        % Plot the updated spectra after removal of the bad epochs / segments
        h = figure;
        pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
        title([subjprefix_plot ' - no bad channels / epochs']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_no_bad_chans_epochs.png']);
        if (par.run_cmd)
            close(h);
        end
    elseif (stage.bad_channels_epochs == ss.fit)
        fprintf('[%s - stage 1 - fit] identifying bad channels and epochs\n', subjprefix_plot);
        EEG = bci2eeglab(ses_data.BCI, 1, par.electrode_file, 1, 1);

        % Remove channels CB1 and CB2
        removed_channels = {'CB1', 'CB2'};
        EEG = pop_select(EEG, 'nochannel', removed_channels);

        % Downsample for speed
        EEG = pop_resample(EEG, 250);
        EEG = eeg_checkset(EEG);

        % Plot the spectra of all channels before any cleaning
        h = figure;
        pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
        title([subjprefix_plot ' - original']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_original.png']);
        if (par.run_cmd)
            close(h);
        end
        
        %% Plot the data to use for visual inspection
        % Band-pass filtering (1-45 Hz) for visual inspection
        if (par.run_manual)
            EEG_filt = EEG;
            [b, a] = butter(2, [1 45] / (EEG_filt.srate / 2)); % define filter parameters 
            EEG_filt.data = filtfilt_boundary(b, a, double(EEG_filt.data)', EEG_filt.event)'; % apply filter
            EEG_filt = eeg_checkset(EEG_filt); % check data for inconsistencies
            pop_eegplot(EEG_filt, 1, 0, 1, [], 'winlength', 20, 'spacing', 50);
        end

        %% Exclude bad epochs
        if (par.run_manual)
            [channels_excluded, trials_excluded, hs] = extract_spectra_outliers(EEG, ...
                [1 1 15; 45 15 45], 3, 0.05, 0.05, 1, par.laplacian_channels, 1);
        else
            [channels_excluded, trials_excluded, hs] = extract_spectra_outliers_auto(EEG, ...
                [1 1 15; 45 15 45], 3, 0.05, 0.05, par.laplacian_channels, 1, 1, 1, 1);
        end
        
        % Export the iteration snapshots and close them
        for i = 1:numel(hs)
            exportgraphics(hs(i), [savepath subjfolder subjprefix '_bad_chanep_iteration_' num2str(i) '.png']);
            close(hs(i));
        end  

        % Remove the bad trials to continue cleaning the spectra
        if ~isempty(trials_excluded)
            epochs_excluded = trials2epochs(EEG, trials_excluded);
            EEG = pop_select(EEG, 'nopoint', epochs_excluded);
            EEG = eeg_checkset(EEG); % check data for inconsistencies
            if isempty(EEG.data)
                warning(['Excluded all data for ' subjprefix_plot]);
                return
            end
        end

        %% Exclude bad channels
        % Use info about bad channels from clean_artifacts
        [~, ~, ~, bad_channels_clean_artifacts] = clean_artifacts(EEG, ...
            'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, ...
            'LineNoiseCriterion', 4, 'Highpass', [0.25 0.75], ...
            'BurstCriterion', 'off', 'WindowCriterion', 'off', ...
            'BurstRejection', 'off');
        bad_channels_clean_artifacts = find(bad_channels_clean_artifacts);
        fprintf('Bad channels as identified by clean_artifacts:\n');
        disp(bad_channels_clean_artifacts');

        % Display bad channels identified by the recursive procedure
        fprintf('Bad channels identified by the recursive procedure:\n');
        disp(channels_excluded');

        % Join all candidates for being a bad channel
        channels_excluded = unique([channels_excluded;
            bad_channels_clean_artifacts]);
        fprintf('All candidates for being a bad channel:\n');
        disp(channels_excluded');

        % plot spectra of the bad channels to reject / confirm them
        nfft = 512;
        noverlap = nfft / 2;
        maxfreq = 50;
        [spec, f] = pwelch(EEG.data', nfft, noverlap, nfft, EEG.srate);
        spec = 10 * log10(spec);
        f_disp = f <= maxfreq;
        mkdir([savepath subjfolder bad_chans_prefix]);
        mkdir([savepath subjfolder good_chans_prefix]);
        is_channel_good = true(EEG.nbchan, 1);
        is_channel_good(channels_excluded) = false;
        for ch = 1:EEG.nbchan
            h = figure; hold on;
            % plot other channels first
            for i = 1:EEG.nbchan
                if (i ~= ch)
                    p = plot(f(f_disp), spec(f_disp, i), 'k');
                    p.Color = [0.8 0.8 0.8];
                end
            end
            % plot the bad channel last so its spectra is always visible
            plot(f(f_disp), spec(f_disp, ch), 'k');
            title([EEG.chanlocs(ch).labels ' (' num2str(ch) ')']);

            good_or_bad = good_chans_prefix;
            if ~is_channel_good(ch)
                good_or_bad = bad_chans_prefix;
            end
            exportgraphics(h, [savepath subjfolder good_or_bad subjprefix '_chan_' num2str(ch) '_' EEG.chanlocs(ch).labels '.png']);
            
            if (par.run_cmd)
                close(h);
            end
        end
        
        if (par.run_manual)
            %%% BREAKPOINT - VISUAL ANALYSIS OF BAD CHANNELS %%%
            channels_corrected = input('Corrected channels:\n');        
        else
            channels_corrected = channels_excluded;
        end       

        % Update the removed_channels field
        removed_channels = [removed_channels {EEG.chanlocs(channels_corrected).labels}];
        disp(removed_channels);

        % remove bad channels
        if ~isempty(channels_corrected)
            EEG = pop_select(EEG, 'nochannel', channels_corrected);
            EEG = eeg_checkset(EEG); % check data for inconsistencies
            
            if isempty(EEG.data)
                warning(['Excluded all data for ' subjprefix_plot]);
                return
            end
        end

        % plot the updated spectra after removal of the bad epochs / segments
        h = figure;
        pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
        title([subjprefix_plot ' - no bad channels / epochs']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_no_bad_chans_epochs.png']);
        if (par.run_cmd)
            close(h);
        end
        
        % helper function for visual analysis on demand
        % pop_eegplot(EEG, 1, 0, 1, [], 'winlength', 20, 'spacing', 50);

        % save info about bad epochs
        preproc_info.bad_trials = trials_excluded;
        preproc_info.bad_channels = removed_channels;
        save([auxpath preproc_info_filename], 'preproc_info');
        fprintf('Preprocessing information saved:\n\n');
        disp(preproc_info);
        disp(preproc_info.bad_trials);
    else
        error('Unexpected stage status');
    end


    %% Stage 2: average reference, filtering in the 1-45 Hz for visual inspection and ICA 
    if (stage.hp_averef == ss.apply)
        fprintf('[%s - stage 2 - apply] applying high-pass filter and average reference\n', subjprefix_plot);
        
        % Average reference
        EEG = pop_reref(EEG, [], 'keepref', 'off');

        % 1 Hz high-pass filter
        [b,a] = butter(2, 1 / (EEG.srate / 2), 'high'); % define filter parameters 
        EEG.data = filtfilt_boundary(b, a, double(EEG.data)', EEG.event)'; % apply filter between boundary events
        EEG = eeg_checkset(EEG); % check data for inconsistencies

        h = figure;
        pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
        title([subjprefix_plot ' - average reference']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_hp_averef.png']);
        if (par.run_cmd)
            close(h);
        end
    elseif (stage.hp_averef == ss.stop)
        % to be completed if needed
    else
        error('Unexpected stage status');
    end


    %% Stage 3: apply ICA
    if (stage.ica == ss.apply)
        fprintf('[%s - stage 3 - apply] copying saved ICA weights\n', subjprefix_plot);
        
        % Copy the ICA weights
        EEG.icaweights = preproc_info.icaweights;
        EEG.icasphere = preproc_info.icasphere;
        EEG.icachansind = preproc_info.icachansind;
        EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);
    elseif (stage.ica == ss.fit)
        fprintf('[%s - stage 3 - fit] applying ICA\n', subjprefix_plot);
        
        % run ICA
        EEG = pop_runica(EEG, 'fastica', 'approach', 'symm');

        % Save the ICA weights
        preproc_info.icaweights = EEG.icaweights;
        preproc_info.icasphere = EEG.icasphere;
        preproc_info.icachansind = EEG.icachansind;
        save([auxpath preproc_info_filename], 'preproc_info');
        fprintf('Preprocessing information was updated with the results of ICA\n');
    elseif (stage.ica == ss.stop)
        % to be completed if needed
    else
        error('Unexpected stage status');
    end

    %% Stage 4: apply ICLabel, check the components and reject artifactual ones
    if (stage.ica_reject == ss.apply)
        fprintf('[%s - stage 4 - apply] rejecting previously identified bad ICA components\n', subjprefix_plot);
        
        % Copy the rejected components
        EEG.reject.gcompreject = preproc_info.gcompreject;

        % Remove the rejected components
        EEG = pop_subcomp(EEG, [], 0, 0);
    elseif (stage.ica_reject == ss.fit)
        fprintf('[%s - stage 4 - fit] identifying bad ICA components\n', subjprefix_plot);
        EEG = iclabel(EEG);
        % Mild reject - only high probability artifacts
        EEG = pop_icflag(EEG, [0 0;           % Brain
                               0.8 1;         % Muscle
                               0.8 1;         % Eye
                               0.8 1;         % Heart
                               0.8 1;         % Line Noise
                               0.8 1;         % Channel Noise
                               0.8 1]);       % Other
        reject_mild = EEG.reject.gcompreject;
        % Hard reject - everything with low probability of being a brain
        % component is rejected
        EEG = pop_icflag(EEG, [0 0.2;         % Brain
                               0.8 1;         % Muscle
                               0.8 1;         % Eye
                               0.8 1;         % Heart
                               0.8 1;         % Line Noise
                               0.8 1;         % Channel Noise
                               0.8 1]);       % Other
        reject_hard = EEG.reject.gcompreject;

        % Get amount of components that explain 95% of the variance
        num_ica_comps = size(EEG.icawinv, 2);
        l = 1; r = num_ica_comps;
        thresh = 95;
        while (l < r - 1)
            m = floor(0.5 * (l + r));
            [~, pvaf] = compvar(EEG.data, {EEG.icasphere, EEG.icaweights}, EEG.icawinv, 1:m);
            if (pvaf < thresh)
                l = m;
            else
                r = m;
            end
        end
        num_comps_95 = r;
        [~, pvaf_95] = compvar(EEG.data, {EEG.icasphere, EEG.icaweights}, ...
                               EEG.icawinv, 1:num_comps_95);
        disp([num_comps_95 pvaf_95]);

        % Use hard reject for the component with high variance and mild for all
        % the others
        reject_joint = reject_mild;
        reject_joint(1:num_comps_95) = reject_hard(1:num_comps_95);
        EEG.reject.gcompreject = reject_joint;
        bad_comps = find(reject_joint');
        good_comps = find(~reject_joint');

        % Plot spectra of bad (to be rejected) and good (to be kept)
        % components
        icaact = eeg_getdatact(EEG, 'component', 1:num_ica_comps);
        
        h = figure; 
        spectopo(icaact, EEG.pnts, EEG.srate, 'freqrange', [1 100], ...
            'plotchans', bad_comps, 'title', [subjprefix_plot ' - bad components']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_bad_ICA_comps.png']);
        if (par.run_cmd)
            close(h);
        end
        
        h = figure; 
        spectopo(icaact, EEG.pnts, EEG.srate, 'freqrange', [1 100], ...
            'plotchans', good_comps, 'title', [subjprefix_plot ' - good components']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_good_ICA_comps.png']);
        if (par.run_cmd)
            close(h);
        end     
        
        % plot spectra of the individual components
        nfft = 512;
        noverlap = nfft / 2;
        maxfreq = 50;
        [spec_ica, f] = pwelch(icaact', nfft, noverlap, nfft, EEG.srate);
        spec_ica = 10 * log10(spec_ica);
        f_disp = f <= maxfreq;
        mkdir([savepath subjfolder bad_comps_prefix]);
        mkdir([savepath subjfolder good_comps_prefix]);
        for ch = 1:num_ica_comps
            h = figure; hold on;
            % plot other channels first
            for i = 1:num_ica_comps
                if (i ~= ch)
                    p = plot(f(f_disp), spec_ica(f_disp, i), 'k');
                    p.Color = [0.8 0.8 0.8];
                end
            end
            % plot the bad channel last so its spectra is always visible
            plot(f(f_disp), spec_ica(f_disp, ch), 'k');
            title(['ICA component ' num2str(ch)]);

            good_or_bad = good_comps_prefix;
            if reject_joint(ch)
                good_or_bad = bad_comps_prefix;
            end
            exportgraphics(h, [savepath subjfolder good_or_bad subjprefix '_ICA_comp_' num2str(ch) '.png']);
            
            if (par.run_cmd)
                close(h);
            end
        end

        if (par.run_manual)
            %%% BREAKPOINT - VISUAL ANALYSIS OF BAD COMPONENTS %%%
%             pop_eegplot(EEG, 1, 0, 1, [], 'winlength', 20, 'spacing', 50);
            pop_viewprops(EEG, 0, bad_comps, {'electrodes', 'off'});
            pop_viewprops(EEG, 0, good_comps, {'electrodes', 'off'});
            fprintf('Pause: reject or accept the components');
            pause();
            % reject_joint(comp) = 0|1; for manual correction
        end

        % Reject only the corrected components
        EEG.reject.gcompreject = reject_joint;
        preproc_info.gcompreject = reject_joint;
        preproc_info.reject_hard = reject_hard;
        preproc_info.reject_mild = reject_mild;
        
        % Remove bad components from the data
        cleanEEG = pop_subcomp(EEG, [], 0, 0);
        
        % Plot the spectra of the result
        h = figure;
        pop_spectopo(cleanEEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
        title([subjprefix_plot ' - without bad ICA components']);
        exportgraphics(h, [savepath subjfolder subjprefix '_spectra_after_ICA.png']);
        if (par.run_cmd)
            close(h);
        end
        
        % Plot the topographies of rejected and kept components
        bad_comps = find(reject_joint');
        pop_topoplot(EEG, 0, bad_comps, [subjprefix_plot ' - bad ICA components'], ...
            0, 0, 'iclabel', 'on');
        h = gcf;
        h.WindowState = 'maximized';
        pause(1);  % wait for the figure to maximize, takes some time if many components are displayed
        exportgraphics(h, [savepath subjfolder subjprefix '_topomaps_bad_ICA_comps.png']);
        if (par.run_cmd)
            close(h);
        end
        
        good_comps = find(~reject_joint');
        pop_topoplot(EEG, 0, good_comps, [subjprefix_plot ' - good ICA components'], ...
            0, 0, 'iclabel', 'on');
        h = gcf;
        h.WindowState = 'maximized';
        pause(1);  % wait for the figure to maximize, takes some time if many components are displayed
        exportgraphics(h, [savepath subjfolder subjprefix '_topomaps_good_ICA_comps.png']);
        if (par.run_cmd)
            close(h);
        end  

        % Save the rejected components
        save([auxpath preproc_info_filename], 'preproc_info');
        fprintf('Saved the rejected components\n');

        % Save the intermediate state
        pop_saveset(EEG, 'filepath', savedata, 'filename', [subjprefix '.set']);
    else
        error('Unexpected stage status');
    end

    %% Stage 5: find the remaining bad epochs based on the peak2peak threshold
    if (stage.peak2peak_reject == ss.fit)
        EEG_epoch = pop_epoch(EEG, {'1'}, [-2.0 8.04]);
        max_amp = max(EEG_epoch.data, [], 2);
        min_amp = min(EEG_epoch.data, [], 2);
        peak2peak = double(squeeze(max_amp - min_amp));
        zpeak2peak = zscore(peak2peak);
        % TODO: correct for already rejected epochs
        zscore_reject = sum(zpeak2peak > 5, 1) > 5;
        absolute_reject = sum(peak2peak > 200, 1);
        bad_epoch_idx = find(zscore_reject | absolute_reject);
    elseif (stage.peak2peak_reject == ss.skip)
        % skipping
    else
        error('Unexpected stage status');
    end
end
     