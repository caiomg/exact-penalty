function [is_eactive, lb_active, ub_active] = ...
        l1_identify_constraints(cmodel, x, lb, ub, epsilon)
    
    tol_bounds = sqrt(eps(1));
    
    epsilon_bounds = min(tol_bounds, epsilon/16);
    
    is_eactive = (abs([cmodel.c]) < epsilon)';
    lb_active = x - lb < epsilon_bounds;
    ub_active = ub - x < epsilon_bounds;
    
end
