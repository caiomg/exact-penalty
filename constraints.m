function [c, ceq] = constraints(c_in, ceq, x)

n_in = length(c_in);
c = zeros(n_in, 1);
for k = 1:n_in
    c(k) = c_in{k}(x);
end

n_eq = length(ceq);
ceq = zeros(n_eq, 1);
for k = 1:n_eq
    c(k) = ceq{k}(x);
end

end