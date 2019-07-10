function [q1, q2, q3, d, Q, R, ind_qr, multipliers, lb_active, ub_active] = ...
        l1_measure_criticality_new(fmodel, cmodel, mu, epsilon, x, lb, ub, threshold)
% L1_MEASURE_CRITICALITY_NEW - 
%   
    tol_g = 1e-10;
    dim = size(x, 1);
    
    [con_eactive, lb_active, ub_active] = l1_identify_constraints(cmodel, ...
                                                      x, lb, ub, epsilon);
    [Q_ext, R_ext, ind_qr, lb_active, ub_active] = ...
        l1_factor_constraint_space(cmodel, con_eactive, lb_active, ...
                                   ub_active);
    p_grad = l1_pseudo_gradient(fmodel.g, mu, cmodel, ind_qr, true);

    try
    n_constraints = numel(ind_qr);
    n_bounds = sum(lb_active) + sum(ub_active);
    q1 = l1_projected_gradient(p_grad, Q_ext, R_ext);
    if q1 > threshold
        [Q, R] = qrdelete_fix(Q_ext, R_ext, n_constraints + 1:n_constraints+ n_bounds);
        q2 = 0;
        q3 = 0;
        d = [];
        multipliers = [];
        % no need to mark bounds
        lb_active = false(dim, 1);
        ub_active = false(dim, 1);
    else
        [degenerate, d, con_confirmed, lb_confirmed, ub_confirmed] ...
            = l1_test_degeneracy(p_grad, Q_ext, R_ext, con_eactive, lb_active, ub_active, mu);
        if degenerate
            multipliers = [];
            q2 = 0;
            q3 = -p_grad'*d;
            current_cols = [con_eactive', lb_active', ub_active'];
            all_confirmed = [con_confirmed', lb_confirmed', ub_confirmed'];
            %all_confirmed = [con_confirmed', false(1, dim), false(1, dim)];
            cols_not_confirmed = current_cols & ~all_confirmed;
            remove_cols = as_col_index(current_cols, cols_not_confirmed);
            [Q, R] = qrdelete_fix(Q_ext, R_ext, remove_cols);
            ind_qr(remove_cols(remove_cols <= numel(ind_qr))) = [];
            lb_active = lb_confirmed;
            ub_active = ub_confirmed;
        else
            q3 = 0;
            [multipliers, lb_multipliers, ub_multipliers, tol_mult] = ...
                l1_estimate_multipliers_new(p_grad, Q_ext, R_ext, ...
                                            lb_active, ub_active);
            [Q_ext, R_ext, lb_active, ub_active] = ...
                l1_remove_bound_constraints(Q_ext, R_ext, con_eactive, ...
                                            lb_active, ub_active, ...
                                            lb_multipliers, ub_multipliers);
            n_bounds = sum(lb_active) + sum(ub_active);
            q1 = l1_projected_gradient(p_grad, Q_ext, R_ext);
            if q1 > threshold
                q2 = 0;
                d = [];
                remove_cols = [n_constraints+1: ...
                               n_constraints+n_bounds];
                [Q, R] = qrdelete_fix(Q_ext, R_ext, remove_cols);
            else
                [q2, d, removed] = l1_choose_search_direction_new( ...
                    Q_ext, R_ext, multipliers, p_grad, mu, x, ...
                        lb, ub, tol_mult, threshold);
                if q2 > q1 ...
                      && 0 < removed && removed <= n_constraints
                        remove_cols = [removed, n_constraints+1: ...
                                       n_constraints+n_bounds];
                        ind_qr = ind_qr([1:removed-1,removed+1:end]);
                else
                    remove_cols = [n_constraints+1: ...
                                   n_constraints+n_bounds];
                end
                [Q, R] = qrdelete_fix(Q_ext, R_ext, remove_cols);
            end
        end
        
    end
    catch myerror
        rethrow(myerror)
    end
    
end
