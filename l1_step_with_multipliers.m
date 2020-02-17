function xm = l1_step_with_multipliers(fmodel_x, cmodel_x, x, ...
                                       multipliers, is_eactive, mu, ...
                                       tr_center, radius, lb, ub)
% L1_STEP_WITH_MULTIPLIERS - 
%   

    dim = numel(x);
    
    g = l1_pseudo_gradient_general(fmodel_x, cmodel_x, mu, [], is_eactive);
    H = l1_pseudo_hessian(fmodel_x, cmodel_x, mu, is_eactive, multipliers);
    
    % Bounds
    lb_shifted = max(lb - x, tr_center - x - radius);
    ub_shifted = min(ub - x, tr_center - x + radius);

    % Orthogonality to nearly-active constraints
    J = [cmodel_x(is_eactive).g]';
    z = zeros(sum(is_eactive), 1);
    
    h = solve_quadratic_problem(H, g, 0, [], [], J, z, lb_shifted, ub_shifted, ...
                                zeros(dim, 1));
    xm = x + h;
    if max(0, xm - ub) + max(0, lb - x) > eps
        'Debug this';
    end
    if norm(xm - tr_center, inf) - radius > eps
        'Debug this';
    end
        
    xm = project_to_bounds(xm, lb, ub);
    
end

