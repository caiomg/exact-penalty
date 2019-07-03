function [Q, R, lb_active, ub_active] = ...
    l1_remove_bound_constraints(Q, R, con_eactive, lb_active, ub_active, lb_multipliers, ub_multipliers)

    tol = 0;

    n_nlconstraints_cols = sum(con_eactive);
    n_lower_cols= sum(lb_active);
    n_upper_cols = sum(ub_active);

    lb_cols_number = n_nlconstraints_cols+1:n_nlconstraints_cols+n_lower_cols;
    ub_cols_number = n_nlconstraints_cols+n_lower_cols+1:n_nlconstraints_cols+n_lower_cols + n_upper_cols;

    lb_active_multiplier = lb_multipliers(lb_active);
    ub_active_multiplier = ub_multipliers(ub_active);
    
    
    lb_remove_col = lb_cols_number(lb_active_multiplier <= tol);
    ub_remove_col = ub_cols_number(ub_active_multiplier <= tol);

    [Q, R] = qrdelete_fix(Q, R, [lb_remove_col, ub_remove_col]);

    lb_active = lb_multipliers > tol;
    ub_active = ub_multipliers > tol;

end
