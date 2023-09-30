function [EEG] = bci2eeglab(BCI, trials_or_task, electrode_file, trial_events, crop_rest)
% EXP2EEGLAB Import experimental data and events to EEGLAB (should be
% initialized behorehand)
% 
% Parameters:
%   BCI - data from the original dataset
%   trials_or_task - indices of trials or the number of task to select
%     if empty, defaults to all trials (1-450)
%     if scalar (task), trials are selected based on BCI.TrialData
%     otherwise used as is (array with indices of trials)
%   electrode_file - file with electrode locations
%   trial_events - whether events with trial IDs should be added (10XXX, 
%     where XXX is the number of the trials)
%   crop_rest - crop repeating segments of the data if identified
%
% Returns:
%   EEG - EEGLAB structure with continuous data (trials are concatenated,
%   breaks in the recording are marked with boundary events)

    % Resolve trials_or_task
    if isempty(trials_or_task)
        % use all trials by default
        trials = 1:numel(BCI.TrialData);
    elseif isscalar(trials_or_task)
        % select trials for a certain task if one number is specified
        trials = find([BCI.TrialData(:).tasknumber] == trials_or_task);
    else
        % interpret as the list of trials to use
        trials = trials_or_task;
    end
    
    % Prepare the data, crop duplicate rest periods if needed
    if crop_rest
        session_data = [];
        n_trials = numel(trials);
        boundary_events = [];
        prev_after_fb = [];
        samples_to_check = 1000;
        for t = trials            
            resultind = BCI.TrialData(t).resultind;
            session_data = cat(2, session_data, ...
                               BCI.data{t}(:, 1:resultind-1));
            rest_after_fb = BCI.data{t}(:, 1:samples_to_check);  % take the data in the beginning of the trial
                           
            % If the beginning of the current trial and the end of the
            % previous trial are not the same, there was a break in the
            % recording
            if ~isempty(prev_after_fb) && ~isequal(rest_after_fb, prev_after_fb)
                boundary_events = [boundary_events t];
            end
            
            prev_after_fb = BCI.data{t}(:, resultind:resultind+samples_to_check-1); % take the data after fb end
        end
    else
        session_data = cat(2, BCI.data{trials});
    end
    
    % Import data
    nbchan = numel(BCI.chaninfo.label);
    EEG = pop_importdata('data', double(session_data), 'srate', BCI.SRATE, 'nbchan', nbchan);
    
    % Load channel locations
    EEG.chanlocs = struct('labels', [BCI.chaninfo(:).label]);
    EEG = pop_chanedit(EEG, 'lookup', electrode_file);
    assert(isequal(BCI.chaninfo.label(:)', {EEG.chanlocs(:).labels}), ...
        "Channel order was corrupted during conversion to EEGLAB");
    
    % Load the information about noisy channels into chanlocs structure
    for ch = 1:EEG.nbchan
        EEG.chanlocs(ch).noisechan = ismember(ch, BCI.chaninfo.noisechan);
    end
    
    % Import events
    events = zeros(3 * numel(trials), 6);
    c = 0;
    offset = 0;
    for t = trials
        if (trial_events)
            % Add 10000 to the trial number to avoid collision with other
            % triggers
            typ = [0, 10000 + t, 1, 2, 3];
            lat = [1, 501, 2001, 4001, BCI.TrialData(t).resultind - 1];
        else
            typ = [0, 1, 2, 3];
            lat = [1, 2001, 4001, BCI.TrialData(t).resultind - 1];
        end
        if ismember(t, boundary_events)
            c = c + 1;
            events(c, :) = [-1, offset, NaN, NaN, NaN, mod(t, 25) == 1];
        end
        for e = 1:numel(typ)
            c = c + 1;
            events(c, :) = [typ(e), offset + lat(e), BCI.TrialData(t).tasknumber, ...
                BCI.TrialData(t).targetnumber, BCI.TrialData(t).result, mod(t, 25) == 1];
        end
        if crop_rest
            offset = offset + BCI.TrialData(t).resultind - 1;
        else
            offset = offset + length(BCI.time{t});
        end
    end
    assert(offset == size(EEG.data, 2), 'Offset should match EEG length');

    EEG = pop_importevent(EEG, 'event', events, 'fields', ...
                          {'type', 'latency', 'task', 'target', 'result', 'new_run'}, ...
                          'append', 'no', 'timeunit', 1 / BCI.SRATE);
    for event_id = find([EEG.event(:).type] == -1)
        EEG.event(event_id).type = 'boundary';
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
                      
    % Import feedback values, trial and meta data
    EEG.etc.positionx = cat(2, BCI.positionx{trials});
    EEG.etc.positiony = cat(2, BCI.positiony{trials});
    EEG.etc.trialdata = BCI.TrialData(trials);
    EEG.etc.metadata = BCI.metadata;
    EEG.etc.chaninfo = BCI.chaninfo;
end