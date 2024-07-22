%%%
% Export connectivity and SNR values for statistical analysis in R
%
% Outputs:
% 1. *output_filename* (provided outside of the script) with the following 
% variables:
% 
%   data - all the results in the long format
%   labels - column names
%   multiverse_labels - labels of different pipelines
%   inv_methods - labels of inverse methods
%   bands - labels of bands
%   roi_methods - labels of ROI aggregation methods (1SVD / 3SVD / AVG-flip)
%   roi_agg_labels - labels of ROI extraction methods (same as above
%       combined with mask / no mask)
%%%

analyzed_sessions = cellfun(@(x) ~isempty(x), preproc_info);

control_reshape = cell(n_multiverse, n_subjects, n_sessions);
for i_multiverse = 1:n_multiverse
    for subject = 1:n_subjects
        for session = 1:n_sessions
            if analyzed_sessions(subject, session) && good_sessions(subject, session)
                control_reshape{i_multiverse, subject, session} = ...
                    [i_multiverse subject session task_accuracy(subject, session) multiverse_data(i_multiverse, :)];
            else
                % mark session with NaN accuracy to be removed later
                control_reshape{i_multiverse, subject, session} = ...
                    [i_multiverse subject session NaN multiverse_data(i_multiverse, :)];
            end
        end
    end
end
metadata = cat(1, control_reshape{:});

% Resting-state
conndata = {
    results_conn.absicoh_within, results_conn.absicoh_across, ...
    results_conn.abslagcoh_within, results_conn.abslagcoh_across, ...
    results_conn.coh_within, results_conn.coh_across, ...
};
snrdata = results_snr.snr_avg;
data = cat(2, metadata, conndata{1}(:), conndata{2}(:), conndata{3}(:), ...
    conndata{4}(:), conndata{5}(:), conndata{6}(:), snrdata(:));
data(isnan(data(:, 4)), :) = [];
labels = {'Pipeline', 'Subject', 'Session', 'Accuracy', ...
          'Inverse', 'Band', 'ROI_Agg', 'Mask', 'ROI_Method', ...
          'ImCoh_Within', 'ImCoh_Across', 'LagCoh_Within', ...
          'LagCoh_Across', 'Coh_Within', 'Coh_Across', 'SNR'};
save([savedata output_filename], ...
    'data', 'labels', 'multiverse_labels', 'inv_methods', 'bands', 'roi_methods', 'roi_agg_labels');