
tries = 12;
n_problems = 134
succeeded = false(n_problems, tries);

for k = 1:tries
    succeeded(:, k) = check_run_single(all_results{k*n_problems});
end
