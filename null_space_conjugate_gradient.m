function s = null_space_conjugate_gradient(fmodel, cmodel, mu, Q, ...
                                           R, ind_qr, x0, d0, radius, ...
                                           lb, ub, s0)
    % NULL_SPACE_CONJUGATE_GRADIENT - 
    %

    dim = size(x0, 1);
    if nargin < 12 || isempty(s0)
        s0 = zeros(dim, 1);
    end

    tol_d = eps(1);
    tol_h = eps(1);
    tol_con = 1e-6;

    s = s0;
    x = project_to_bounds(x0 + s, lb, ub);    
    [Q, R, bl_active, bu_active] = detect_and_include_active_bounds(Q, ...
                                                      R, x, d0, lb, ...
                                                      ub, tol_con);
    r_cols = size(R, 2);
    N = Q(:, r_cols+1:end);
    N(bl_active|bu_active, :) = zeros(sum(bl_active | bu_active), ...
                                      dim - r_cols);
    d = N*(N'*d0);

    pred_s = predict_descent(fmodel, cmodel, s, mu, []);

    for iter = 1:dim
        
        s_prev = s;
        d(bl_active | bu_active) = zeros(sum(bl_active | bu_active), 1);
        if norm(d) < tol_d
            break
        end

        fmodel_d = shift_model(fmodel, s);
        for k = 1:length(cmodel)
            cmodel_d(k) = shift_model(cmodel(k), s);
        end
        x = project_to_bounds(x0 + s, lb, ub);
        [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, d);
        tmax_tr = tr_radius_breakpoint(d, radius, s);
        
        tmax = min(tmax_bounds, tmax_tr);
        [t, brpoints_crossed] = line_search_cg(fmodel_d, cmodel_d, mu, d, tmax);
        pred_d = predict_descent(fmodel, cmodel, s + t*d, mu, []);
        if pred_d > pred_s
            s = s + t*d;
        else
            break
        end
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

            
            pgradient = l1_pseudo_gradient_new(fmodel, cmodel, mu, s, ind_qr);
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
                mid_interval = 0.5*interval_length;
                mid_point = s_prev + (bpoint + mid_interval)*d;
                
                H = H + (interval_length/t)*l1_hessian(fmodel, cmodel, ...
                                                       mu, mid_point);
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

    end
end
