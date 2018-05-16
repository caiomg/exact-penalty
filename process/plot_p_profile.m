k = 1
%%
any_solved = sum([cobyla_solved, l1_solved(:, k)], 2) > 0;
computing_time = [l1_evals(:, k), cobyla_evals];
[rho, tau] = dm_performance_profile(computing_time);
figure(2)
semilogx(tau, rho);
figure(1)
perf(computing_time, 1);
k = k+1;
% figure(3)
% plot(tau, rho)