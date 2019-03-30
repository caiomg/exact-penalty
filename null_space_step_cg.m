function s = null_space_step_cg(fmodel, cmodel, mu, x0, ind_qr, Q, R, ...
                                radius, lb, ub)

                            
    tol_d = eps(1);
    tol_h = eps(1);
    tol_con = 1e-6;
    


    pgradient = l1_pseudo_gradient_new(fmodel, cmodel, mu, [], ind_qr);
    [Q, R, bl_active, bu_active] = ...
        detect_and_include_active_bounds(Q, R, x0, -pgradient, lb, ub, tol_con);

    active_bounds = bl_active | bu_active;
    [dim, r_cols] = size(R);
    N = Q(:, r_cols+1:end);
    N(active_bounds, :) = zeros(sum(active_bounds), dim - r_cols);
    g = (N*(N'*pgradient));
    d = -g;
    
    s = zeros(dim, 1);

    fmodel_d = fmodel;
    cmodel_d = cmodel;
    
    for iter = 1:dim

        d(bl_active | bu_active) = zeros(sum(bl_active | bu_active), 1);
        if norm(d) < tol_d
            break
        end
        
        x = project_to_bounds(x0 + s, lb, ub);
        [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, d);
        tmax_tr = tr_radius_breakpoint(d, radius, s);
        
        tmax = min(tmax_bounds, tmax_tr);
        [t, brpoints_crossed] = line_search_cg(fmodel_d, cmodel_d, mu, d, tmax);
        s = s + t*d;
        if tmax_bounds ~= 0 && (t == 0 || t == tmax_tr)
            break
        elseif tmax_bounds == 0 || t == tmax
            % Repeat all those projections
            % Restart method
            bl_active_new = t_lb == t;
            bu_active_new = t_ub == t;
            s(bl_active_new) = lb(bl_active_new) - x0(bl_active_new);
            s(bu_active_new) = ub(bu_active_new) - x0(bu_active_new);
            
            
            [Q, R, bl_included, bu_included] = ...
                     detect_and_include_active_bounds(Q, R, x, d, lb, ub, tol_con);
            bl_active = bl_active | bl_included;
            bu_active = bu_active | bu_included;

            fmodel_d = shift_model(fmodel, s);
            for k = 1:length(cmodel)
                cmodel_d(k) = shift_model(cmodel(k), s);
            end
            
            pgradient = l1_pseudo_gradient_new(fmodel_d, cmodel_d, mu, [], ind_qr);
            [dim, r_cols] = size(R);
            N = Q(:, r_cols+1:end);
            N(bl_active | bu_active, :) = ...
                zeros(sum(bl_active | bu_active), dim - r_cols);
            g = (N*(N'*pgradient));
            d = -g;
        else
            % Stopped at a local minimum
            % Compute conjugate directions
            
            % Compute average Hessian
            bpoint_prev = 0;
            H = zeros(dim);
            for k = 1:length(brpoints_crossed)
                bpoint = brpoints_crossed(k);
                interval_length = (bpoint - bpoint_prev);
                mid = 0.5*interval_length;
                
                H = H + (interval_length/t)*l1_hessian(fmodel_d, cmodel_d, ...
                                                       mu, mid*d);
                bpoint_prev = bpoint;
            end
            
            % Conjugate direction
            pgradient_new = l1_pseudo_gradient_new(fmodel, cmodel, mu, s, ind_qr);
            g_new = N*(N'*pgradient_new);
            g_new(bl_active | bu_active) = ...
                zeros(sum(bl_active | bu_active), 1);
            dHd = d'*H*d;
            if norm(dHd) < tol_h
                break
            end
            beta = (g_new'*H*d)/dHd;
            d = -g_new + beta*d;
        end
        s(bl_active) = lb(bl_active) - x0(bl_active);
        s(bu_active) = ub(bu_active) - x0(bu_active);

        fmodel_d = shift_model(fmodel, s);
        for k = 1:length(cmodel)
            cmodel_d(k) = shift_model(cmodel(k), s);
        end
    end
    
    
end
