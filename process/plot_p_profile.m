k = 1
%%
any_solved = sum([cobyla_solved, l1_solved(:, k+1)], 2) > 0;
[rho, tau] = dm_performance_profile([cobyla_evals, l1_evals(:, k+1)], any_solved);
semilogx(tau, rho);
k = k+1;