function [con_multipliers, lb_multipliers, ub_multipliers, tol] = ...
        l1_estimate_multipliers_with_bounds(p_grad, Q, R, lb_active, ub_active)
% L1_ESTIMATE_MULTIPLIERS_NO_BOUNDS - 
%   

    dim = size(p_grad, 1);
    H = R'*R;
    g = (R'*Q'*p_grad);
    
    f = @(x) quadratic(H, g, 0, x);
    
    n_lower_active = sum(lb_active);
    n_upper_active = sum(ub_active);
    n_bounds_active = n_lower_active + n_upper_active;
    n_cons_total = size(H, 2);
    n_cons_not_bounds = n_cons_total - n_bounds_active;

    bl_mult = [-inf(n_cons_not_bounds, 1);
               zeros(n_bounds_active, 1)];
    
    
    fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                   'Algorithm', 'active-set', ...
                                   'SpecifyObjectiveGradient', true, ...
                                   'OptimalityTolerance', 1e-6);

    [all_multipliers, ~, exitflag, ~, fmincon_lambda] = ...
        fmincon(f, zeros(size(g)), [], [], [], [], bl_mult, [], [], fmincon_options);
    
    con_multipliers = all_multipliers(1:n_cons_not_bounds);
    lb_multipliers = zeros(dim, 1);
    ub_multipliers = zeros(dim, 1);
    
    lb_multipliers(lb_active) = all_multipliers(n_cons_not_bounds+1: ...
                                     n_cons_not_bounds+ n_lower_active);
    ub_multipliers(ub_active) = all_multipliers(n_cons_not_bounds+ n_lower_active+1:end);
    
    assert(norm(all_multipliers(fmincon_lambda.lower > 0), inf) == 0);

    tol = sqrt(eps);

                                     
    
end
