function [passed] = test_sensor2voxel()
    data = [1 1 1 1 1; 2 2 2 2 2; 3 3 3 3 3];
    data_ep = cat(3, data, data + 10, data + 20, data + 30);
    A_inv_1D = [0 1 0; 0 0 1; 1 0 0]'; % rotate
    A_inv_3D = cat(3, [0 1 0; 0 0 1; 1 0 0]', [0 0 1; 1 0 0; 0 1 0]', eye(3));

    data_inv = circshift(data, -1, 1)';
    data_ep_inv = reshape(circshift(data_ep, -1, 1), 3, [])';
    data_ep_inv_smask = data_ep_inv;
    data_ep_inv_smask(:, 1) = 0;
    data_inv_3D = cat(3, circshift(data, -1, 1)', circshift(data, 1, 1)', data');
    data_ep_inv_3D = cat(3, ...
        reshape(circshift(data_ep, -1, 1), 3, [])', ...
        reshape(circshift(data_ep, 1, 1), 3, [])', ...
        reshape(data_ep,  3, [])');

    test_cases = struct( ...
        'data', {data, data_ep, data_ep, data_ep, data, data_ep}, ...
        'smask', {1:3, [1 3], 1:3, 1:3, 1:3, 1:3}, ...
        'A_inv', {A_inv_1D, A_inv_1D([1 3], :), A_inv_1D, A_inv_1D, A_inv_3D, A_inv_3D}, ...
        'vmask', {1:3, 1:3, [1 3], 1:3, 1:3, 1:3}, ...
        'expected_result', {data_inv, data_ep_inv_smask, data_ep_inv(:, [1 3]), data_ep_inv, data_inv_3D, data_ep_inv_3D}, ...
        'desc', {'cnt_1D', 'ep_1D_smask', 'ep_1D_vmask', 'ep_1D', 'cnt_3D', 'ep_3D'});
    num_test_cases = numel(test_cases);

    passed = false(num_test_cases);
    for i = 1:num_test_cases
        tc = test_cases(i);
        fprintf('%s_%s...', mfilename(), tc.desc);
        voxel_data = sensor2voxel(tc.data, tc.smask, tc.A_inv, tc.vmask);
        assert(isequal(voxel_data, tc.expected_result));       
        fprintf('OK\n');
        passed(i) = 1;
    end
end