function [cfg] = BCI_MI_init_workspace(toolboxes, python_path)
%BCI_MI_INIT_WORKSPACE Configure all paths that are necessary
%
% Returns:
%   cfg - structure with the paths to assets, data (raw and preprocessed),
%   results (intermediate and final), and toolboxes

    toolbox.base = normpath([pwd() '/../toolboxes/']);

    asset.base = normpath([pwd() '/../assets/']);

    data.base = normpath([pwd() '/../data/']);
    data.raw = [data.base 'raw/'];
    
    preproc.base = [data.base 'preproc/'];
    preproc.aux = [data.base 'aux/'];
    preproc.task1 = [preproc.base 'task1/'];

    results.base = normpath([pwd() '/../results/']);
    results.misc = [results.base 'miscellaneous/'];
    results.fooof = [results.base 'fooof/'];
    results.fooof_lap = [results.fooof 'laplace/'];
    results.fooof_roi = [results.fooof 'multiverse/'];
    results.tex = [results.base 'tex/'];

    folders_to_create = {preproc.base, preproc.aux, ...
        preproc.task1, results.base, results.misc, ...
        results.fooof, results.fooof_lap, results.fooof_roi, results.tex};
    for i = 1:numel(folders_to_create)
        if ~exist(folders_to_create{i}, 'dir')
           mkdir(folders_to_create{i});
        end
    end
    
    % Add relevant code folders to the path
    p = strsplit(mfilename('fullpath'), '/');
    script_folder = strjoin(p(1:end-1), '/');
    addpath(script_folder);
    folders_to_add = {'connectivity', 'csp', 'inverse', 'laplace', ...
        'multiverse', 'performance', 'precomputed', 'preprocessing', ...
        'sanity', 'snr', 'tests', 'util'};
    for i = 1:numel(folders_to_add)
        addpath([script_folder '/' folders_to_add{i}]);
    end
    
    n_toolboxes = numel(toolboxes);
    for tb = 1:n_toolboxes
        fprintf('Loading toolbox %s...\n', toolboxes{tb});
        switch (toolboxes{tb})
            case 'bbci'
                toolbox.bbci = [toolbox.base 'bbci_public/'];
                toolbox.btb = [toolbox.base 'BTB/BTB.mat'];
                addpath(genpath(toolbox.bbci));
            case 'eeglab'
                toolbox.eeglab = [toolbox.base 'eeglab2021.0/'];
                addpath(toolbox.eeglab);
                eeglab;
            case 'fastica'
                toolbox.fastica = [toolbox.base 'FastICA_25/'];
                addpath(toolbox.fastica);
            case 'fooof'
                toolbox.fooof = [toolbox.base 'fooof_mat/fooof_mat/'];
                addpath(toolbox.fooof);
            case 'haufe'
                toolbox.haufe = [toolbox.base 'haufe/'];
                addpath(toolbox.haufe);
            case 'python'
                toolbox.python = python_path;
    
                % Executing Python code in a separate process protects from crashing
                % the whole MATLAB if Python crashes, but implied a 2 GB limit on the
                % size of shared data (??)
                pe = pyenv;
                if strcmp(pe.Executable, toolbox.python) && pe.Status == 0
                    terminate(pe);
                end

                pyenv('Version', toolbox.python, 'ExecutionMode', 'OutOfProcess');
            case 'tprod'
                % NOTE: tprod tests passed only when everything was in double precision
                toolbox.tprod = [toolbox.base 'tprod/'];
                addpath(genpath(toolbox.tprod));

                tprod([1 2 3], [1 -1], [1; 2; 3], [-1 2]); % test call
            otherwise
                error(['Unknown toolbox ', toolboxes{tb}]);
        end
        fprintf('OK\n');
    end   
    
    %% Pack all the info into one structure
    cfg.asset = asset;
    cfg.data = data;
    cfg.preproc = preproc;
    cfg.results = results;
    cfg.toolbox = toolbox;
end

function [np] = normpath(p)
% Quick solution to normalize relative paths
    d = dir(p);
    np = [d(1).folder '/'];
end
