function p = l1_function(f, phi, mu, x, ind_eactive)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    ind_eactive = [];
end

n_constraints = length(phi);
sum_phi = 0;
for n = 1:n_constraints
    if isempty(find(ind_eactive == n, 1))
        sum_phi = sum_phi + max(phi{n}(x), 0);
    end
end
p = f(x) + mu*sum_phi;

end

