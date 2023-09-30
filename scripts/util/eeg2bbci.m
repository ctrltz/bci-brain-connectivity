function [epo] = eeg2bbci(EEG, labels)
% EEG2BBCI Extract epochs from EEG structure into BBCI format structure  
%
% Parameters:
%   EEG - epoched EEG data
%   labels - class for each epoch
%
% Returns:
%   epo - structure that is suitable for the BBCI toolbox

    epo.fs = EEG.srate;
    epo.clab = {EEG.chanlocs(:).labels};
    epo.x = permute(double(EEG.data), [2 1 3]);
    epo.x = squeeze(epo.x);
    
    if EEG.trials > 1
        y = unique(labels);
        epo.y = double(y' == labels);
    end
end
