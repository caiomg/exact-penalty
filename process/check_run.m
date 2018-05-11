
tries = length(all_results);
n_problems = 134
succeeded = false(n_problems, tries);
all_results2 = {};
for k = 1:tries
    succeeded(:, k) = check_run_single(all_results{k});
    all_results2{end+1} = all_results{k};
end


b = sum(succeeded, 2);
bi = [-max(0, min(1, b));
      zeros(tries, 1)];
c = ones(tries, 1);
intcon = 1:tries;
A = [-succeeded;
     -eye(tries)];

chosen = intlinprog(c, intcon, A, bi)