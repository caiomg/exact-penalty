function p = l1_function(f, phi, mu, x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_constraints = length(phi);
sum_phi = 0;
for n = 1:n_constraints
    sum_phi = sum_phi + min(phi{n}(x), 0);
end
p = f(x) - sum_phi/mu;

end

