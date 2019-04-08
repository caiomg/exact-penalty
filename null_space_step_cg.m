function s = null_space_step_cg(fmodel, cmodel, mu, x0, ind_qr, Q, R, ...
                                radius, lb, ub)


    tol_con = 1e-6;
    
    pgradient = l1_pseudo_gradient_new(fmodel, cmodel, mu, [], ind_qr);
    [Q, R, bl_active, bu_active] = ...
        detect_and_include_active_bounds(Q, R, x0, -pgradient, lb, ub, tol_con);

    d = -pgradient;
    
    s = null_space_conjugate_gradient(fmodel, cmodel, mu, Q, R, ...
                                      ind_qr, x0, d0, radius, lb, ...
                                      ub, bl_active, bu_active);

end
