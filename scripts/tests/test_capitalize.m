function [passed] = test_capitalize()

    test_cases = struct(...
        'name', {'left', 'right'}, ...
        'expected_result', {'Left', 'Right'}, ...
        'desc', {'left', 'right'});
    num_test_cases = numel(test_cases);

    passed = false(num_test_cases, 1);
    for i = 1:num_test_cases
        tc = test_cases(i);
        fprintf('%s_%s...', mfilename(), tc.desc);
        capitalized = capitalize(tc.name);
        assert(strcmp(capitalized, tc.expected_result));       
        fprintf('OK\n');
        passed(i) = 1;
    end
end

