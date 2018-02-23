function [p, fvalues] = l1_function(f, phi, mu, x, ind_eactive)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    ind_eactive = [];
end

n_functions = 1 + length(phi);
fvalues = zeros(n_functions, 1);
value_used = false(n_functions, 1);
fvalues(1) = f(x);
for n = 1:n_functions-1
    fvalues(n+1, 1) = phi{n}(x);
    if isempty(find(ind_eactive == n, 1))
        value_used(n+1) = true;
    end
end
sum_phi = sum(max(fvalues(value_used), 0));
p = fvalues(1) + mu*sum_phi;

end

