%%%
% Merge files with preprocessing info
%
% Outputs:
% 1. BCI_MI_preproc_info.mat with the following variables:
%   preproc_info {subjects x sessions} - preprocessing info for each
%       subject and session if the data was present
%   good_sessions [subjects x sessions] - manually identified good sessions
%%%

auxpath = cfg.preproc.aux;

n_subjects = 62;
n_sessions = 11;
preproc_info = cell(n_subjects, n_sessions);
good_sessions = zeros(n_subjects, n_sessions);
for subject = 1:n_subjects
    for session = 1:n_sessions
        subjprefix = ['S' num2str(subject) '_Session_', num2str(session)];
        preproc_info_filename = [subjprefix '_preproc_info.mat'];
        if exist([auxpath preproc_info_filename], 'file')   
            tmp = load([auxpath preproc_info_filename]);
            preproc_info{subject, session} = tmp.preproc_info;
            preproc_info{subject, session}.bad_trials = ...
                sort(tmp.preproc_info.bad_trials);
            good_sessions(subject, session) = 1;
        end
    end
    fprintf('%d ', subject);
end
fprintf('Done\n');

% manually identified bad sessions (subject / session)
% some recordings contain spikes without a clear pattern so that they
% could be removed manually -> spikes ruin ICA, leak in all frequencies
bad_sessions = [16 1; 31 1; 53 4];
for i = 1:size(bad_sessions, 1)
    subject = bad_sessions(i, 1);
    session = bad_sessions(i, 2);
    good_sessions(subject, session) = 0;
end

save([auxpath 'BCI_MI_preproc_info.mat'], 'preproc_info', 'good_sessions');
