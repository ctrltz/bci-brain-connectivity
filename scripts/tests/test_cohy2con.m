function [passed] = test_cohy2con()
    cohy_test = [1 + 0i, 0 + 1i, -1 + 0i, 0 - 1i, (1 + 1i) / sqrt(2), (-1 + 1i) / sqrt(2), ...
        (-1 - 1i) / sqrt(2), (1 - 1i) / sqrt(2), (sqrt(3) - 1i) / 4, (-1 + sqrt(3)*1i) / 4];
    expected_conn = {
        [0, 1,    0,  -1,    1/sqrt(2), 1/sqrt(2), -1/sqrt(2), -1/sqrt(2), -0.25,      sqrt(3)/4], ...  icoh
        [0, 1,    0,  1,     1/sqrt(2), 1/sqrt(2), 1/sqrt(2),  1/sqrt(2),  0.25,       sqrt(3)/4], ...  absicoh
        [1, 1,    1,  1,     1,         1,         1,          1,          0.5,        0.5      ], ...  coh
        [0, 1,    0,  -1,    1,         1,         -1,         -1,         -0.2773501, 0.4472136], ...  lagcoh
        [0, 1,    0,  1,     1,         1,         1,          1,          0.2773501,  0.4472136], ...  abslagcoh
        [0, pi/2, 0,  -pi/2, pi/4,      -pi/4,     pi/4,       -pi/4,      -pi/6,      -pi/3    ], ...  angle
        [0, pi/2, 0,  pi/2,  pi/4,      pi/4,      pi/4,       pi/4,       pi/6,       pi/3     ], ...  absangle
        [0, pi/2, pi, -pi/2, pi/4,      3*pi/4,    -3*pi/4,    -pi/4,      -pi/6,      2*pi/3   ], ...  angle2
        [0, pi/2, pi, pi/2,  pi/4,      3*pi/4,    3*pi/4,     pi/4,       pi/6,       2*pi/3   ], ...  absangle2
    };
    
    test_cases = struct(...
        'measure', {'icoh', 'icoh', 'coh', 'lagcoh', 'lagcoh', 'angle', 'angle', 'angle2', 'angle2'}, ...
        'use_abs', {0, 1, 0, 0, 1, 0, 1, 0, 1}, ...
        'expected_result', expected_conn, ...
        'desc', {'icoh', 'absicoh', 'coh', 'lagcoh', 'abslagcoh', 'angle', 'absangle', 'angle2', 'absangle2'});
    num_test_cases = numel(test_cases);
    
    passed = false(num_test_cases, 1);
    for i = 1:num_test_cases
        tc = test_cases(i);
        fprintf('%s_%s...', mfilename(), tc.desc);
        conn = cohy2con(cohy_test, tc.measure, tc.use_abs);
        assert(max(abs(conn - tc.expected_result)) < 1e-7);       
        fprintf('OK\n');
        passed(i) = 1;
    end
end