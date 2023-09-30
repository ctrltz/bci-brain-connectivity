function [passed] = test_get_voxels_roi()
    sa_test.cortex75K.in_HO = [1 1 1 1 2 2 2 2 3 3 3 3]';
    sa_test.voxels_5K_cort = [1 3 5 7 9 11]';
    test_voxel_mask = [1 1 1 0 0 1]';
    
    test_cases = struct(...
        'roi_ind', {1 2 3 1 2 3}, ...
        'expected_result', {[1; 2], 3, 6, [1; 2], [3; 4], [5; 6]}, ...
        'mask', {test_voxel_mask, test_voxel_mask, test_voxel_mask, [], [], []}, ...
        'desc', {'1_mask', '2_mask', '3_mask', ...
                 '1_nomask', '2_nomask', '3_nomask'});
    num_test_cases = numel(test_cases);
    
    passed = false(num_test_cases, 1);
    for i = 1:num_test_cases
        tc = test_cases(i);
        fprintf('%s_%s...', mfilename(), tc.desc);
        voxels_roi = get_voxels_roi(sa_test, tc.roi_ind, tc.mask);
        assert(isequal(voxels_roi, tc.expected_result));       
        fprintf('OK\n');
        passed(i) = 1;
    end
end