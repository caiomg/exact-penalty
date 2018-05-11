k = 1
%%
any_solved = sum([cobyla_solved, l1_solved(:, k)], 2) > 0;
computing_time = [cobyla_evals, l1_evals(:, k)];
[rho, tau] = dm_performance_profile(computing_time);
figure(2)
semilogx(tau, rho);
figure(1)
perf(computing_time, 1);
k = k+1;
figure(3)
plot(tau, rho)