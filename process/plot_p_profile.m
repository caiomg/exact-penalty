k = 1
%%
computing_time = [my_evals_cheti, al_evals(:, :)];
[rho, tau] = dm_performance_profile(computing_time);
% figure(2)
semilogx(tau, rho, 'LineWidth', 1.0);
ylim([0, 1]);
legend('$\ell_1$-penalty', 'AL-CSearch',  'AL-BOBYQA', 'AL-NMead', 'Location', 'southeast')
ylabel('$\rho_s(\alpha)$ (Fraction of problems solved by solver)');
xlabel('$\alpha$ (Performance ratio: function evaluations)')
filename = '/home/caio/pdoc/exact-penalty-paper/submission-4/current/figures/performance_al.tikz';
matlab2tikz(filename, 'parseStrings', false, 'width', '0.9\columnwidth')
% figure(1)
% perf(computing_time, 1);
% k = k+1;
% figure(3)
% plot(tau, rho)

%%
computing_time = [my_evals_sampaio, sampaio_evals(:, :)];
[rho, tau] = dm_performance_profile(computing_time);
% figure(2)
semilogx(tau, rho, 'LineWidth', 1.0);
ylim([0, 1]);
legend('$\ell_1$-penalty', 'TF Subbasis', 'TF Min. Fro norm', 'TF Min 2-norm', 'TF Regression', 'Location', 'southeast')
ylabel('$\rho_s(\alpha)$ (Fraction of problems solved by solver)');
xlabel('$\alpha$ (Performance ratio: function evaluations)')
filename = '/home/caio/pdoc/exact-penalty-paper/submission-4/current/figures/performance_tf.tikz';
matlab2tikz(filename, 'parseStrings', false, 'width', '0.9\columnwidth')

%% Sampaio 2
computing_time = [my_evals_sampaio_penalty_known, sampaio_evals(:, :)];
[rho, tau] = dm_performance_profile(computing_time);
% figure(2)
semilogx(tau, rho, 'LineWidth', 1.0);
ylim([0, 1]);
legend('$\ell_1$-penalty', 'TF Subbasis', 'TF Min. Fro norm', 'TF Min 2-norm', 'TF Regression', 'Location', 'southeast')
ylabel('$\rho_s(\alpha)$ (Fraction of problems solved by solver)');
xlabel('$\alpha$ (Performance ratio: function evaluations)')
filename = '/home/caio/pdoc/exact-penalty-paper/submission-4/current/figures/performance_tf_known.tikz';
matlab2tikz(filename, 'parseStrings', false, 'width', '0.9\columnwidth')
