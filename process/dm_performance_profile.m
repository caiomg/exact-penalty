function [rho, tau] = dm_performance_profile(computing_time)

[n_problems, n_solvers] = size(computing_time);


best_time = min(computing_time, [], 2);
success = ~isnan(best_time);
n_success = sum(success);

% Removing failures from all algorithms
computing_time = computing_time(success, :);
best_time = best_time(success);

r = zeros(size(computing_time));
for k = 1:n_success
    r(k, :) = computing_time(k, :)/best_time(k);
end

tau = unique(sort(r(:)));
tau = tau(~isnan(tau)); % removing rM
n_tau = length(tau);
rho = zeros(n_tau, n_solvers);

for k = 1:n_tau
    t = tau(k);
    for s = 1:n_solvers
        rho(k, s) = sum(r(:, s) <= t)/n_problems;
    end
end

end
