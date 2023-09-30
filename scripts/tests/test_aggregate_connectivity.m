function [passed] = test_aggregate_connectivity()
    conn_matrix1 = repmat([
        11 12 13 14;
        12 22 23 24;
        13 23 33 34;
        14 24 34 44
	], 3, 3);
    conn_matrix1 = reshape(conn_matrix1, 1, 12, 12);

    conn_matrix2 = [
        11 12 13 14;
        12 22 23 24;
        13 23 33 34;
        14 24 34 44
    ];
    expected_result = conn_matrix2;
    conn_matrix2 = reshape(conn_matrix2, 1, 4, 4);
    
    band_freqs = false(20, 1);
    band_freqs(5:8) = true;
    conn_matrix3 = zeros(20, 12, 12);
    conn_matrix3(band_freqs, :, :) = repmat(conn_matrix1, sum(band_freqs), 1, 1);
    
    conn_matrix4 = zeros(20, 4, 4);
    conn_matrix4(band_freqs, :, :) = repmat(conn_matrix2, sum(band_freqs), 1, 1);
    
    test_cases = struct(...
        'data', {{conn_matrix2}, {conn_matrix1}, {conn_matrix1, conn_matrix1}, ...
                 {conn_matrix4}, {conn_matrix3}, {conn_matrix3, conn_matrix3}}, ... 
        'n_comps', {1 3 3 1 3 3}, ...
        'n_rois', {4 4 4 4 4 4}, ...
        'band_freqs', {[], [], [], band_freqs, band_freqs, band_freqs}, ...
        'expected_result', {expected_result, expected_result, {expected_result, expected_result}, ...
                            expected_result, expected_result, {expected_result, expected_result}}, ...
        'desc', {'nb_1comp', 'nb_3comps', 'nb_3comps_multiple', ...
                 'bb_1comp', 'bb_3comps', 'bb_3comps_multiple'});
    num_test_cases = numel(test_cases);

    passed = false(num_test_cases, 1);
    for i = 1:num_test_cases
        tc = test_cases(i);
        fprintf('%s_%s...', mfilename(), tc.desc);
        n_data = numel(tc.data);
        if n_data == 1
            conn_agg = aggregate_connectivity(tc.data, tc.n_comps, tc.n_rois, tc.band_freqs);
            assert(isequal(conn_agg, tc.expected_result));
        else
            [conn_agg1, conn_agg2] = aggregate_connectivity(tc.data, tc.n_comps, tc.n_rois, tc.band_freqs);
            assert(isequal(conn_agg1, tc.expected_result{1}));
            assert(isequal(conn_agg2, tc.expected_result{2}));
        end
        fprintf('OK\n');
        passed(i) = 1;
    end
end