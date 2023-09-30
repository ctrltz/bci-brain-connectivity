function [EEG, EEG_narrow] = prepare_data(EEG, all_chanlocs, freq_band, downsample, n_mirror_pnts)
% PREPARE_DATA Preprocess the data before actual calculations (interpolate
% missing channels, downsample, and filter in a narrow band)
%
% Parameters:
%   EEG - data to be preprocessed
%   all_chanlocs - structure with locations of all channels
%   freq_band - frequency band for narrow-band data
%   downsample - new sampling frequency
%   n_mirror_pnts - number of samples to use for mirroring edge effects of
%       filtering
%
% Returns:
%   EEG - preprocessed data (broad-band)
%   EEG_narrow - preprocessed data (narrow-band)

    % Interpolate bad channels that were removed from the dataset    
    [~, bad_elec] = setdiff({all_chanlocs.labels}, {EEG.chanlocs.labels});
    if numel(bad_elec)
        EEG = pop_interp(EEG, all_chanlocs, 'spherical');
    end
    
    % Down-sample the data to speed up some of the computations
    EEG = pop_resample(EEG, downsample);
    
    % Apply a narrow-band filter
    EEG_narrow = EEG;
    [b, a] = butter(2, 2 * freq_band / EEG_narrow.srate);
    EEG_narrow.data = filtfilt_boundary(b, a, double(EEG_narrow.data'), EEG.event, n_mirror_pnts)';
end
