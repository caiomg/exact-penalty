function [measure, d, is_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel_x, ...
                                                     cmodel_x, x, ...
                                                     mu, epsilon, ...
                                                     lb, ub, center, radius)
% L1_CRITICALITY_MEASURE_AND_DESCENT_DIRECTION - 
%   

    global as_succeeded
    global ip_succeeded
    global as_lower
    global total_measure_difference
    if isempty(as_succeeded)
        as_succeeded = 0;
        ip_succeeded = 0;
        as_lower = 0;
        total_measure_difference = 0;
    end
    
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
    
        d1 = solve_linear_problem(f, Aineq, bineq, [], [], plb, pub);
        [measure, d] = correct_measure_computation(pg, G, mu, dlb, dub, d1);

end
