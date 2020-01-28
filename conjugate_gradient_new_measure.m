function x = conjugate_gradient_new_measure(fmodel, cmodel, x0, mu, epsilon, ...
                                            radius, lb, ub, d0)
    % NULL_SPACE_CONJUG
    %

    dim = size(x0, 1);

    tol_d = eps(1);
    tol_h = eps(1);
    tol_con = 1e-6;
    
    fmodel_d = fmodel;
    cmodel_d = cmodel;

    x = x0;
    if nargin < 8
            [~, d] = l1_criticality_measure_and_descent_direction(fmodel_d, ...
                                                              cmodel_d, ...
                                                              x, mu, ...
                                                              epsilon, ...
                                                              lb, ...
                                                              ub, x0, radius);
    else
        d = d0;
    end

    if ~isempty(find(x > ub | x < lb, 1))
        warning('cmg:out_of_bounds', 'Point already out of bounds');
    end

    for iter = 1:dim
        
        if norm(d) < tol_d
            break
        end
        
        [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, d);
        try
        tmax_tr = infinity_tr_radius_breakpoint(d, radius, x - x0);
        catch teste
            rethrow(teste);
        end
    
        tmax = min(tmax_bounds, tmax_tr);
        [t, brpoints_crossed] = line_search_cg(fmodel_d, cmodel_d, mu, d, tmax);
    
        lower_bound_hits = t_lb == t;
        upper_bound_hits = t_ub == t;
        x = x + t*d;
        x(lower_bound_hits) = lb(lower_bound_hits);
        x(upper_bound_hits) = ub(upper_bound_hits);
        
        s = x - x0;
        if ~isempty(find(abs(s) >= radius, 1))
            break
        end
        fmodel_d = shift_model(fmodel, s);
        for k = 1:length(cmodel)
            cmodel_d(k) = shift_model(cmodel(k), s);
        end
        try
        [sigma, g_neg] = l1_criticality_measure_and_descent_direction(fmodel_d, ...
                                                          cmodel_d, ...
                                                          x, mu, ...
                                                          epsilon, ...
                                                          lb, ub, x0, radius);
        catch teste
            rethrow(teste);
        end
        
        if ~isempty(find(upper_bound_hits | lower_bound_hits, 1))
            d = g_neg; % Restarting
        else
            % Compute average Hessian
            bpoint_prev = 0;
            H = zeros(dim);
            for k = 1:length(brpoints_crossed)
                bpoint = brpoints_crossed(k);
                interval_length = (bpoint - bpoint_prev);
                mid_interval = 0.5*interval_length;
                mid_point = x + (bpoint + mid_interval)*d;
                
                H = H + (interval_length/t)*l1_hessian(fmodel, cmodel, ...
                                                       mu, mid_point);
                bpoint_prev = bpoint;
            end
            
            % Conjugate direction
            dHd = d'*H*d;
            if norm(dHd) < tol_h
                break
            end
            beta = ((g_neg'*H*d)/dHd);
            d = g_neg + beta*d;
            
        end
    end
end
