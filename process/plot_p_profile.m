k = 1
%%
computing_time = [my_evals, al_evals(:, :)];
[rho, tau] = dm_performance_profile(computing_time);
% figure(2)
semilogx(tau, rho, 'LineWidth', 1.0);
ylim([0, 1]);
legend('$\ell_1$-penalty', 'AL-CSearch',  'AL-BOBYQA', 'AL-NMead', 'Location', 'southeast')
ylabel('$\rho_s(\alpha)$ (Fraction of problems solved by solver)');
xlabel('$\alpha$ (Performance ratio: function evaluations normalized by best)')
matlab2tikz('~/docd/paper_sources_new/figures/performance3.tikz', 'parseStrings', false, 'width', '0.85\columnwidth')
% figure(1)
% perf(computing_time, 1);
% k = k+1;
% figure(3)
% plot(tau, rho)