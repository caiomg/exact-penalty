
tries = 12;
n_problems = 134
succeeded = false(n_problems, tries);

for k = 1:tries
    succeeded(:, k) = check_run_single(all_results{k*n_problems});
end


b = sum(succeeded, 2);
bi = [-max(0, min(1, b));
      zeros(tries, 1)];
c = ones(tries, 1);
intcon = 1:tries;
A = [-succeeded;
     -eye(tries)];

chosen = intlinprog(c, intcon, A, bi)