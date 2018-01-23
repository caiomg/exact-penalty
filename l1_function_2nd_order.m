function p = l1_function_2nd_order(f, phi, mu, x, ind_eactive, multipliers, ind_multipliers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    ind_eactive = [];
end

n_constraints = length(phi);
sum_c = 0;
sum_m = 0;
for n = 1:n_constraints
    sum_c = sum_c + max(phi{n}(x), 0);
    if ~isempty(find(ind_multipliers == n, 1))
        sum_m = sum_m + multipliers(ind_multipliers == n)*phi{n}(x);
    end
end
p = f(x) + mu*sum_c + sum_m;

end

