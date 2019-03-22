function v = l1_range_step_complete(fmodel, cmodel, Q, R, mu, h, ind_qr, ...
                                    radius, x0, lb, ub)
    
    [dim, r_cols] = size(R);
    lb_active = false(dim, 1);
    ub_active = false(dim, 1);
    tol_radius = 0.05*radius + max(0, norm(h) - radius);

    radius_icb = radius;    

    I = eye(dim);
    
    c_vals = update_constraint_information(cmodel, ind_qr, h);
    
    x = x0 + h;
    At = R'*Q;
    Q_range = Q(:, 1:r_cols);
    while true
        Z = I(:, ~(lb_active | ub_active));
        Q_range_bounds = orth(Z*(Z'*Q_range));
        if size(Q_range_bounds, 2) == 0
            break
        end
        % Set TR problem
        At_red = At*Q_range_bounds;
        H = At_red'*At_red;
        g = At_red'*c_vals;
        % Compute direction
        try
        s_red = ms_step(H, g, radius_icb);
        catch thiserror
            rethrow(thiserror)
        end
        % Direction in full space
        s_full = Q_range_bounds*s_red;
        s_full(lb_active) = zeros(sum(lb_active), 1); % Shouldn't
        s_full(ub_active) = zeros(sum(ub_active), 1); % change much

        % Check bounds
        [lb_active_s, ub_active_s] = active_bounds(x, s_full, lb, ub);

        lb_active_new = lb_active_s & ~lb_active;
        ub_active_new = ub_active_s & ~ub_active;
        lb_active = lb_active | lb_active_new;
        ub_active = ub_active | ub_active_new;
        
        if norm((x + s_full) - x0) > radius + tol_radius
            radius_icb = 0.5*radius_icb;
        elseif isempty(find(lb_active_new | ub_active_new, 1))
            [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, s_full);
            tmax_tr = tr_radius_breakpoint(s_full, radius, x - x0);

            sHs = s_red'*H*s_red;
            gs = g'*s_red;
            step_change = 0.5*sHs + gs;
            if tmax_tr >= 1 && tmax_bounds >= 1
                x = x + s_full;
                break
            elseif tmax_tr <= tmax_bounds ...
                          && 0.5*sHs*tmax_tr^2 + gs*tmax_tr < 0.5*step_change
                x = x + tmax_tr*s_full;
                break
            elseif tmax_bounds <= tmax_tr ...
                    && 0.5*sHs*tmax_bounds^2 + gs*tmax_bounds < 0.5*step_change
                x = x + tmax_bounds*s_full;

                lb_active_new = t_lb == tmax_bounds;
                ub_active_new = t_ub == tmax_bounds;
                lb_active = lb_active | lb_active_new;
                ub_active = ub_active | ub_active_new;
            else
                radius_icb = 0.5*radius_icb;
            end
        end
    end
    v = x - (x0 + h);
    
end
