function [] = BCI_MI_preprocessing_apply(subject, session, prepinfo, par)
% BCI_MI_PREPROCESSING_APPLY Apply preprocessing using precomputed
% information about bad trials, channels, and ICA components
%
% Parameters:
%   subject - subject ID (1-62)
%   session - session ID (1-11)
%   prepinfo - structure that contains preprocessing information
%   par - additional parameters (paths, flags that control what to save)

    %% Unpack par
    datapath = par.data.raw;
    savedata = par.preproc.base;
    savedata_raw = par.savedata_raw_eeglab;

    %% Load data
    filename = ['S' num2str(subject) '_Session_' num2str(session) '.mat'];
    if ~exist([datapath filename], 'file')
        warning([filename ' - file not found, skipping']);
        return
    end
    
    fprintf('Loading data from %s\n', filename);    
    ses_data = load([datapath filename]);
    BCI = ses_data.BCI;
    fprintf('Done\n');
    
    %% Plot spectra of the raw data and save if needed
    EEG = bci2eeglab(BCI, 1, par.electrode_file, 1, 1);
    if (par.save_raw)
        savefile = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        pop_saveset(EEG, 'filepath', [savedata_raw savefile]);
    end
    
    h = figure;
    subplot(1, 2, 1);
    pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
    title('raw data');
    
    %% Add information about bad trials and extract EEG 
    % Data from task 1 was analyzed
    task = 1;

    % Extract good trials for a given task
    task_trials = [BCI.TrialData(:).tasknumber] == task;
    task_trialdata = BCI.TrialData(task_trials);

    good_trials = ~ismember(1:numel(task_trialdata), prepinfo.bad_trials);
    good_task_trials = [task_trialdata(good_trials).trialnumber];
    EEG = bci2eeglab(BCI, good_task_trials, par.electrode_file, 1, 1);
    
    % Check that all trials belong to the specified task
    assert(all([EEG.etc.trialdata(:).tasknumber] == task));
    
    % Remove previously identified bad channels
    if ~isempty(prepinfo.bad_channels)
        EEG = pop_select(EEG, 'nochannel', prepinfo.bad_channels);
        EEG = eeg_checkset(EEG); % check data for inconsistencies
    end

    % Re-reference to average reference
    EEG = pop_reref(EEG, [], 'keepref', 'off');

    % Copy the ICA weights
    EEG.icaweights = prepinfo.icaweights;
    EEG.icasphere = prepinfo.icasphere;
    EEG.icachansind = prepinfo.icachansind;
    EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);

    % Copy the rejected components
    EEG.reject.gcompreject = prepinfo.gcompreject;

    % Remove the rejected components
    EEG = pop_subcomp(EEG, [], 0, 0);

    % Remove the DC offset
    is_boundary = strcmp({EEG.event(:).type}, 'boundary');
    boundary_lat = [EEG.event(is_boundary).latency];
    cont_segments = [0 boundary_lat; boundary_lat EEG.pnts];
    n_segments = size(cont_segments, 2);

    fprintf('Removing DC offset from %d continuous segments...', n_segments);
    for i = 1:n_segments
        s_idx = round(cont_segments(1, i)) + 1;
        e_idx = round(cont_segments(2, i));
        
        EEG.data(:, s_idx:e_idx) = double(EEG.data(:, s_idx:e_idx)) - mean(double(EEG.data(:, s_idx:e_idx)), 2);
    end
    fprintf('Done\n');
    
    % Plot the spectra
    subplot(1, 2, 2);
    pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freqrange', [1 100], 'electrodes', 'off');
    ylim([-inf 30]); % focus on the alpha peak, ylim get stretched due to DC offset sometimes
    title(['preprocessed data - task ', num2str(task)]);
    
    % Save the preprocessed data
    if (par.save_preproc)
        savefile = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        pop_saveset(EEG, 'filepath', [savedata 'task' num2str(task) '/' savefile]);
    end
    
    h.WindowState = 'maximized';
    pause(1);   % wait for window to maximize
    sgtitle(['Subject ' num2str(subject) ' - Session ' num2str(session)]);
    exportgraphics(h, [par.savepath 'apply/S' num2str(subject, '%02d') '_Session_' num2str(session, '%02d') '.png']);
    if (par.run_cmd)
        close(h);
    end
end
