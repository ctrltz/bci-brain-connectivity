function [EEG_condition] = extract_condition_EEG(EEG, EEG_narrow, condition, events, windows)
% EXTRACT_CONDITION_EEG Extract EEG data for different conditions (rest,
% preparation, feedback) and classes (right, left)
%
% Parameters:
%   EEG - broadband data
%   EEG_narrow - narrowband data
%   condition - a cell array with 3 elements
%       1) trial class: 1 (right) and/or 2 (left)
%       2) task period: 'rest', 'prep', or 'fbend'
%       3) frequency band: 'bb' (broad) or 'nb' (narrow)
%   events - triggers to use as onset
%   windows - time windows of task periods
%
% Returns:
%   EEG_condition - EEG structure with epoched data according to the
%   specified condition

    %% Extract condition data
    trial_class = condition{1};
    task = condition{2};
    band = condition{3};
    
    % Check the values
    if ~ismember(band, {'bb', 'nb'})
        error('bad band');
    end
    if ~ismember(task, {'rest', 'prep', 'fbend'})
        error('bad task');
    end
    if any(~ismember(trial_class, [1, 2]))
        error('bad trial class');
    end
        
    % Select the data source
    if strcmp(band, 'bb')
        EEG_to_use = EEG;
    else
        EEG_to_use = EEG_narrow;
    end
        
    % Select time period and reference event
    event = {events.target};
    switch (task)
        case 'rest'
            time_period = windows.rest;
        case 'prep'
            time_period = windows.prep;
        case 'fbend'
            time_period = windows.fbend;
            event = {events.fbend};
    end

    % Select trials to pick
    ses_trialdata = EEG.etc.trialdata;
    trials_to_pick = ismember([ses_trialdata(:).targetnumber], trial_class);
        
    disp({trial_class, task, band, event{1}, time_period});
    EEG_epoched = pop_epoch(EEG_to_use, event, time_period);
    EEG_epoched = pop_select(EEG_epoched, 'trial', find(trials_to_pick));
    EEG_condition = eeg_checkset(EEG_epoched);
end
