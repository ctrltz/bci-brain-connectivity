%%%
% Extract metadata and performance from all sessions
%
% Outputs:
%
% 1. BCI_MI_performance.mat with the following variables:
%
%   perfs {subjects x sessions} - information about each session (accuracy,
%       trial length, number of good / bad trials)
%   metadata {subjects x sessions} - metadata from the dataset
%
% 2. BCI_MI_task1_accuracy.mat with the variable:
%
%   task1_accuracy [subjects x sessions] - BCI accuracy
%
% 3. BCI_MI_performance_metadata_long.csv - performance and metadata
% exported to CSV format for analysis in R
%
%%%

%% Main loop - extracting performance
perfs = cell(n_subjects, n_sessions);
metadata = cell(n_subjects, n_sessions);

parfor (subject = 1:n_subjects, num_workers)
% for subject = 1:5
    for session = 1:n_sessions
%     for session = 1
        perf = [];
        perf.Subject = subject;
        perf.Session = session;

        %% Load the original data
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.mat'];
        fprintf('Processing Subject %i, Session %i...\n', subject, session);

        fprintf('Loading data...');
        if ~exist([cfg.data.raw filename], 'file')
            warning("Skipping subject %d session %d - file not found\n", subject, session);
            continue
        end
        ses_data = load([cfg.data.raw filename]);
        fprintf('Done\n');

        BCI = ses_data.BCI;

        % Count the number of bad trials
        bt = preproc_info{subject, session}.bad_trials;
        perf.BadTrials = numel(bt);

        %% Get performance based on all trials
        task_trials = [BCI.TrialData(:).tasknumber] == task;
        task_fbtime = [BCI.TrialData(task_trials).triallength];
        task_result = [BCI.TrialData(task_trials).result];
        task_forcedresult = [BCI.TrialData(task_trials).forcedresult];
        completed_trials = ~isnan(task_result);

        perf.TotalTrials = sum(task_trials);
        perf.MeanTrialLength = mean(task_fbtime);
        perf.ValidTrials = sum(completed_trials);
        perf.MeanValidTrialLength = mean(task_fbtime(completed_trials));
        perf.Accuracy = mean(task_result(completed_trials));
        perf.AccuracyForced = mean(task_forcedresult);

        % Count NaN result as a bad one
        task_result(isnan(task_result)) = 0;
        perf.AccuracyNanIsBad = mean(task_result);

        %% Load EEG without any bad trials
        filename = ['S' num2str(subject) '_Session_' num2str(session) '.set'];
        EEG = pop_loadset('filepath', datapath, 'filename', filename);
        EEG_epoched = pop_epoch(EEG, {'1'}, [-2 2]);

        % Check that the trial information matches the expectations
        assert(numel(EEG.etc.trialdata) == EEG_epoched.trials);  % same number of trials
        assert(all([EEG.etc.trialdata(:).tasknumber] == task));  % all trials are from task 1
        
        %% Get performance based on good trials
        good_task_result = [EEG.etc.trialdata(:).result];
        good_task_fbtime = [EEG.etc.trialdata(:).triallength];
        good_completed_trials = ~isnan(good_task_result);    

        perf.GoodTotalTrials = numel(EEG.etc.trialdata);
        perf.GoodTrialLength = mean(good_task_fbtime);
        perf.GoodValidTrials = sum(good_completed_trials);
        perf.GoodValidTrialLength = mean(good_task_fbtime(good_completed_trials));
        perf.AccuracyGoodTrials = mean(good_task_result(good_completed_trials));
        
        % Report class balance for trials w/o artifacts
        good_completed_targets = [EEG.etc.trialdata(:).targetnumber];
        for i_class = 1:n_classes
            perf.(['ValidTrials' capitalize(classes{i_class})]) = ...
                sum(good_completed_targets == classes_to_analyze(i_class));
        end

        % Count NaN result as a bad one
        good_task_result(isnan(good_task_result)) = 0;
        perf.AccuracyGoodNanIsBad = mean(good_task_result);

        %% Get performance in other tasks for comparison
        for other_task = 1:3
            task_trials = [BCI.TrialData(:).tasknumber] == other_task;
            task_result = [BCI.TrialData(task_trials).result];
            completed_trials = ~isnan(task_result);
            perf = setfield(perf, ['AccuracyTask' num2str(other_task)], ...
                mean(task_result(completed_trials)));
        end

        %% Save the data
        perfs{subject, session} = perf;
        metadata{subject, session} = EEG.etc.metadata;
    end
end

%% Save the results in the original format
save([savedata 'BCI_MI_performance.mat'], 'perfs', 'metadata');

%% Save accuracy values in a separate file
% Create accuracy array
task_accuracy = NaN(n_subjects, n_sessions);
for subject = 1:n_subjects
    for session = 1:n_sessions
        if ~isempty(perfs{subject, session})
            task_accuracy(subject, session) = perfs{subject, session}.Accuracy;
        end
    end
end

save([savedata 'BCI_MI_task_accuracy.mat'], 'task_accuracy');

%% Save the results in the long format
perfs_long = perfs(:);
perfs_long = perfs_long(cellfun(@(x) ~isempty(x), perfs_long));
perfs_long = cell2mat(perfs_long);
T1 = struct2table(perfs_long);

assert(all(([perfs_long(:).GoodTotalTrials] + ...
    [perfs_long(:).BadTrials]) == [perfs_long(:).TotalTrials]));

metadata_long = metadata(:);
metadata_long = metadata_long(cellfun(@(x) ~isempty(x), metadata_long));
metadata_long = cell2mat(metadata_long);

% Fix NaNs in string columns since writetable processes them weirdly
cols = {'handedness', 'instrument', 'athlete', 'handsport', 'hobby', 'gender'};
for col = cols
    nan_cells = cellfun(@(x) isnan(x), {metadata_long(:).(col{1})});
    [metadata_long(nan_cells).(col{1})] = deal({'NaN'});
end

T2 = struct2table(metadata_long);

T = [T1 T2];
writetable(T, [savedata 'BCI_MI_performance_metadata_long.csv']);