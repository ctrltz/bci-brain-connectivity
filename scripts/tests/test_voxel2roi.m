function [passed] = test_voxel2roi()
    test_voxel_data = [1 2 1; 2 1 2; 3 3 3]';
    test_signflip = [1 1 -1]';

    test_cases = struct(...
        'agg_method', {'svd', 'svd', 'avg', 'avg-flip'}, ...
        'n_comps', {1, 3, 1, 1}, ...
        'signflip', {[], [], [], test_signflip}, ...
        'expected_roi_data', {[-0.408248; 0.816496; -0.408248], [], [2; 2; 2], [0; 0; 0]}, ...
        'expected_ev', {1, [1; 0; 0], 0, 0}, ...
        'expected_w', {[0.707106; -0.707106; 0], [], [], []}, ...
        'desc', {'1svd', '3svd', '1avg', '1avg_flip'});
    num_test_cases = numel(test_cases);

    passed = false(num_test_cases, 1);
    for i = 1:num_test_cases
        tc = test_cases(i);
        fprintf('%s_%s...', mfilename(), tc.desc);
        [roi_data, ev, w] = voxel2roi(test_voxel_data, tc.agg_method, tc.n_comps, tc.signflip);
        
        if ~isempty(tc.expected_roi_data)
            assert(max(abs(roi_data - tc.expected_roi_data), [], 'all') < 1e-5);
        else
            assert(size(roi_data, 2) == tc.n_comps);
            assert(size(roi_data, 1) == size(test_voxel_data, 1));
        end
        assert(max(abs(ev - tc.expected_ev), [], 'all') < 1e-5);
        if ~isempty(tc.expected_w)
            assert(max(abs(w - tc.expected_w), [], 'all') < 1e-5);
        elseif strcmp(tc.agg_method, 'svd')
            assert(size(w, 2) == tc.n_comps);
            assert(size(w, 2) == size(test_voxel_data, 2));
        end
        fprintf('OK\n');
        passed(i) = 1;
    end
end

