function xv = l1_feasibility_correction(cmodel, ind_eactive, tr_center, xh, radius, lb, ub)
% L1_FEASIBILITY_CORRECTION - 
%   
    
    h = xh - tr_center;
    [c, A] = constraint_values_and_gradients_at_point(cmodel, ind_eactive, h);
    
    H = A*A';
    g = A*c;
    f = @(v) quadratic(H, g, 0, v);
    
    v_lb = max(-radius - h, lb - tr_center - h);
    v_ub = min(radius - h, ub - tr_center - h);
    
    
    v0 = 0.5*(v_lb + v_ub);

    fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                   'Algorithm', 'active-set', ...
                                   'SpecifyObjectiveGradient', true, ...
                                   'OptimalityTolerance', 1e-6);

    [v, ~, exitflag] = fmincon(f, v0, [], [], [], [], v_lb, v_ub, ...
                               [], fmincon_options);
    xv = xh + v;
    
end
