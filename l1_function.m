function [p, fvalues] = l1_function(f, phi, con_lb, con_ub, mu, x)

    assert(nargin == 6);

    n_constraints = length(phi);
    f_val = f(x);
    con_vals = zeros(n_constraints, 1);
    sum_violations = 0;
    for n = 1:n_constraints
        con_vals(n) = phi{n}(x);
    end
    viol_lb = max(0, con_lb - con_vals);
    viol_ub = max(0, con_vals - con_ub);

    p = f_val+ mu*sum(viol_lb + viol_ub);
    fvalues = [f_val;
               con_vals];

end

