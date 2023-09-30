function [] = BCI_MI_run_tests()
    p = strsplit(mfilename('fullpath'), '/');
    folder = strjoin(p(1:end-1), '/');
    test_files = dir([folder '/test_*.m']);
    
    for i = 1:numel(test_files)
        run(test_files(i).name);
    end
end