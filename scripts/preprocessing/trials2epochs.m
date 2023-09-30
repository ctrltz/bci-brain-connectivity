function [epochs] = trials2epochs(EEG, trials)
% TRIALS2EPOCHS Get start and end latencies for given trials
%
% Parameters:
%   EEG - data in EEGLAB 
%   trials - indices of trials
%
% Returns:
%   epochs [trials x 2] - start and end latency for each trial

    epochs = zeros(numel(trials), 2);
    trial_start = EEG.event(ismember({EEG.event(:).type}, {'0'}));
    trial_end = EEG.event(ismember({EEG.event(:).type}, {'3'}));
    for i = 1:numel(trials)
        tr = trials(i);
        slat = floor(trial_start(tr).latency);
        elat = ceil(trial_end(tr).latency);
        epochs(i, :) = [slat, elat];
    end
end

