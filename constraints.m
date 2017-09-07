function [c, ceq] = constraints(c_in, ceq, x, factor)

if nargin < 4
    factor = 1;
end

n_in = length(c_in);
c = zeros(n_in, 1);
for k = 1:n_in
    c(k) = factor*c_in{k}(x);
end

n_eq = length(ceq);
ceq = zeros(n_eq, 1);
for k = 1:n_eq
    c(k) = factor*ceq{k}(x);
end

end