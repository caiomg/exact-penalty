function s = null_space_step_cg(fmodel, cmodel, mu, x0, ind_qr, Q, R, ...
                                radius, lb, ub, multipliers)

                            
    tol_d = eps(1);
    tol_h = eps(1);
    suff_radius = 0.9*radius;
    [dim, r_cols] = size(R);
    N = Q(:, r_cols+1:end);
    
    pgradient = l1_pseudo_gradient_new(fmodel, cmodel, mu);

    s = zeros(dim, 1);
    g = (N*(N'*pgradient));
    d = -g;
    bl_active = false(dim, 1);
    bu_active = false(dim, 1);
    
    for checks = 1:dim
        [bl_active_x, bu_active_x] = active_bounds(x0, d, lb, ub);
        bl_active_new = bl_active_x & ~bl_active;
        bu_active_new = bu_active_x & ~bu_active;

        if isempty(find(bl_active_new | bu_active_new, 1))
            break
        else
            [Q, R, bl_included, bu_included] = ...
                include_bounds_gradients(Q, R, bl_active_new, ...
                                         bu_active_new);
            r_cols = size(R, 2);
            N = Q(:, r_cols+1:end);
            bl_active = bl_active | bl_active_x;
            bu_active = bu_active | bu_active_x;

            g = N*(N'*pgradient);
            d = -g;
            d(bl_active | bu_active) = ...
                zeros(sum(bl_active | bu_active), 1);
        end
    end
    fmodel_d = fmodel;
    cmodel_d = cmodel;
    
    for iter = 1:dim
        
        if norm(d) < tol_d
            break
        end

        [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x0 + s, lb, ub, d);
        tmax_tr = tr_radius_breakpoint(d, radius, s);
        
        tmax = min(tmax_bounds, tmax_tr);
        [t, brpoints_crossed] = line_search_cg(fmodel_d, cmodel_d, mu, d, tmax);
        s = s + t*d;
        if t == 0 || t == tmax_tr
            break
        elseif t == tmax
            % Repeat all those projections
            % Restart method
            bl_active_new = t_lb == t;
            bu_active_new = t_ub == t;
            
            [bl_active_x, bu_active_x] = active_bounds(x0 + s, d, lb, ub);
            bl_active_new_2 = bl_active_x & ~bl_active;
            bu_active_new_2 = bu_active_x & ~bu_active;
    
            if ~isempty(find((bl_active_new ~= bl_active_new_2) ...
                    | (bu_active_new ~= bu_active_new_2), 1))
                warning('cmg:runtime_error', 'This needs debugging');
            end
            
            for checks = 1:sum(~(bl_active | bu_active))
                if isempty(find(bl_active_new | bu_active_new, 1))
                    if checks == 1
                        warning('cmg:runtime_error', 'This needs debugging');
                    end                        
                    break
                else
                    [Q, R, bl_included, bu_included] = ...
                        include_bounds_gradients(Q, R, bl_active_new, ...
                                                 bu_active_new);
                    r_cols = size(R, 2);
                    N = Q(:, r_cols+1:end);
                    bl_active = bl_active | bl_included;
                    bu_active = bu_active | bu_included;
                    pgradient = l1_pseudo_gradient_new(fmodel, cmodel, mu, s);
                    g = N*(N'*pgradient);
                    d = -g;
                    d(bl_active | bu_active) = ...
                        zeros(sum(bl_active | bu_active), 1);
                end
            end
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
            pgradient_new = l1_pseudo_gradient_new(fmodel, cmodel, mu, s);
            g_new = N*(N'*pgradient_new);
            g_new(bl_active | bu_active) = ...
                zeros(sum(bl_active | bu_active), 1);
            dHd = d'*H*d;
            if norm(dHd) < tol_h
                break
            end
            beta = (g_new'*H*d)/dHd;
            d = -g_new + beta*d;
            d(bl_active | bu_active) = ...
                zeros(sum(bl_active | bu_active), 1);

        end
        x = x0 + s;
        if max([max(0, lb(bl_active) - x(bl_active));
                max(0, x(bu_active) - ub(bu_active))]) > sqrt(eps)
            warning('cmg:bounds_error', 'This needs debugging');
        end

        fmodel_d = shift_model(fmodel, s);
        for k = 1:length(cmodel)
            cmodel_d(k) = shift_model(cmodel(k), s);
        end
    end
    
    
end
