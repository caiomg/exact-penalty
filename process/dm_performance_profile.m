function [rho, tau] = dm_performance_profile(computing_time, success)

[n_problems, n_solvers] = size(computing_time);

if nargin < 2 || isempty(success)
    success = true(n_problems, 1);
end

r = zeros(size(computing_time)) + max(max(computing_time)) + 1;


for k = 1:n_problems
    if success(k)
        r(k, :) = computing_time(k, :)/min(computing_time(k, :));
    end
end

tau = unique(sort(r(:)));
tau = tau(1:end-1); % removing rM
n_tau = length(tau);
rho = zeros(n_tau, n_solvers);

for k = 1:n_tau
    t = tau(k);
    for s = 1:n_solvers
        rho(k, s) = sum(r(:, s) <= t)/n_problems;
    end
end

end
