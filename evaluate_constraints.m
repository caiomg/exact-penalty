function cx = evaluate_constraints(constraints, x, con_lb, con_ub)

    if nargin == 2 || (isempty(con_lb) && isempty(con_ub))
        con_ub = zeros(size(x));
    end
    
    n_constraints = length(constraints);
    ind = 0;
    for n = n_constraints:-1:1
        [phi_c, phi_g, phi_H] = constraints{n}(x);
        if ~isempty(con_ub) && isfinite(con_ub(n))
            ind = ind + 1;
            cx(ind, 1).c = phi_c - con_ub(n);
            cx(ind, 1).g = phi_g;
            cx(ind, 1).H = phi_H;
        end
        if ~isempty(con_lb) && isfinite(con_lb(n))
            ind = ind + 1;
            cx(ind, 1).c = con_lb(n) - phi_c;
            cx(ind, 1).g = -phi_g;
            cx(ind, 1).H = -phi_H;
        end
    end

end
