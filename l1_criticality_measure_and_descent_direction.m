function [sigma, d, is_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel_x, ...
                                                     cmodel_x, x, ...
                                                     mu, epsilon, ...
                                                     lb, ub, center, radius)
% L1_CRITICALITY_MEASURE_AND_DESCENT_DIRECTION - 
%   

    if nargin < 8
        center = x;
        radius = inf;
    end

    dim = size(x, 1);
    is_eactive = l1_identify_constraints(cmodel_x, x, lb, ub, epsilon);

    s0 = zeros(dim, 1);
    pg = l1_pseudo_gradient_general(fmodel_x, cmodel_x, mu, s0, is_eactive);
    
    n_eactive = sum(is_eactive);
    
    f = [pg;
         mu*ones(n_eactive, 1)];
    
    G = [cmodel_x(is_eactive).g];
    I = eye(n_eactive);

    Aineq = [G', -I];
    bineq = zeros(n_eactive, 1);
    
    dlb = max(-1, max(lb - x, (center - x) - radius));
    dub = min(1, min(ub - x, (center - x) + radius));
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
    if exitflag < 0
        'Debug this';
    end
    if sigma_neg > 0
        sigma_neg = 0;
    end
    if sigma_neg > 0 %OLD UNUSED CODE
        if ~isempty(dy) && (isempty(Aineq) || isempty(find(Aineq*(-dy) > bineq, 1)))
            dy = -dy;
            sigma_neg = -sigma_neg;
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
                error('cmg:criticality_error', ...
                      'Could not compute criticality measure');
            end
            if sigma_neg > 0 ...
                 && norm(max(0, Aineq*dy - bineq)) >= norm(max(0, Aineq*(-dy) - bineq))
                dy = project_to_bounds(-dy, dlb, dub);
                sigma_neg = f'*dy;
            end
        end
    end

    if isempty(dy)
        1;
    end
    d = dy(1:dim);
    sigma = -sigma_neg;
    
end
