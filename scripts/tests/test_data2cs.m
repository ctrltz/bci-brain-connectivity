function [passed] = test_data2cs()
% This test is instable but passes after rerunning snapshot_test_data2cs.m
    rng('default');

    p = strsplit(mfilename('fullpath'), '/');
    folder = strjoin(p(1:end-1), '/');
    load([folder '/snapshot_test_data2cs.mat'], 'data', 'fres', 'cs_bb', 'cs_nb');
    
    passed = false(2, 1);

    % bb avg
    fprintf('%s_bb_avg...', mfilename());
    cs = data2cs_bb(data, fres, 1);
    assert(isequal(cs, cs_bb));
    fprintf('OK\n');
    passed(1) = 1;

    % nb avg
    fprintf('%s_nb_avg...', mfilename());
    cs = data2cs_nb(data, 1);
    assert(isequal(cs, cs_nb));
    fprintf('OK\n');
    passed(2) = 1;
end