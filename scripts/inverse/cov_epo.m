function [C] = cov_epo(data_epo, smask)
% COV_EPO Calculate covariance matrix of the epoched data (epochs are 
% concatenated before calculating the covariance)
%
% Parameters: 
%   data_epo [channels x samples x epochs] - EEG data
%   smask - mask to use only a subset of channels
%
% Returns:
%   C - covariance matrix
    
    C = cov(data_epo(smask, :)');
end