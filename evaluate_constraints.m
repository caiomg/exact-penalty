function cx = evaluate_constraints(constraints, x)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
n_constraints = length(constraints);
for n = n_constraints:-1:1
    [phi_c, phi_g, phi_H] = constraints{n}(x);
    cx(n, 1).c = phi_c;
    cx(n, 1).g = phi_g;
    cx(n, 1).H = phi_H;
end

end

