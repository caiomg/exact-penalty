function [sigma, d, ind_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel, cmodel, ...
                                                     mu, x, epsilon, lb, ub)
% L1_CRITICALITY_MEASURE_AND_DESCENT_DIRECTION - 
%   

    dim = size(x, 1);
    ind_eactive = l1_identify_constraints(cmodel, x, lb, ub, epsilon);

    s0 = zeros(dim, 1);
    pg = l1_pseudo_gradient_general(fmodel, cmodel, mu, s0, ind_eactive);
    
    n_eactive = sum(ind_eactive);
    
    f = [pg;
         mu*ones(n_eactive, 1)];
    
    G = [cmodel(ind_eactive).g];
    I = eye(n_eactive);
    Z = zeros(dim, n_eactive);

    Aineq = [G', -I];
    bineq = zeros(n_eactive, 1);
    
    dlb = max(-1, lb - x);
    dub = min(1, ub - x);
    ylb = zeros(n_eactive, 1);
    yub = inf(n_eactive, 1);
    
    plb = [dlb;
           ylb];
    pub = [dub;
           yub];
    
    linprog_problem.solver = 'linprog';
    linprog_problem.options = optimoptions('linprog', 'Display', ...
                                           'off', 'Algorithm', 'dual-simplex');
    linprog_problem.f = f;
    linprog_problem.Aineq = Aineq;
    linprog_problem.bineq = bineq;
    linprog_problem.lb = plb;
    linprog_problem.ub = pub;
    [dy, sigma_neg, exitflag, output] = linprog(linprog_problem);
    if sigma_neg > 0
        if isempty(find(Aineq*(-dy) > bineq, 1))
            dy = -dy;
            sigma_neg = - sigma_neg;
        else
            Aineq = [Aineq;
                    f'];
            bineq = [bineq;
                    0];
            linprog_problem.Aineq = Aineq;
            linprog_problem.bineq = bineq;
            linprog_problem.options = optimoptions('linprog', ...
                                                   'Display', 'off', ...
                                                   'ConstraintTolerance', 1e-7, ...
                                                   'Algorithm', 'dual-simplex');
            [dy, sigma_neg, exitflag, output] = linprog(linprog_problem);
            if exitflag < 0
                warning('bla');
            end
            if sigma_neg > 0
                warning('bla');
            end
        end
    end

    d = dy(1:dim);
    sigma = -sigma_neg;
    
end
