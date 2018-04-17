
tries = length(all_results);
n_problems = length(all_problems)
succeeded = false(n_problems, tries);

for k = 1:tries
    succeeded(:, k) = check_run_single(all_results{k});
end
