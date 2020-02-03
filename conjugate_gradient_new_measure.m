function x = conjugate_gradient_new_measure(fmodel, cmodel, x0, mu, epsilon, ...
                                            radius, lb, ub)
    % NULL_SPACE_CONJUG
    %

    dim = size(x0, 1);
    n_constraints = numel(cmodel);
    max_iters = dim + n_constraints;

    tol_d = eps(1);
    tol_h = eps(1);
    tol_con = 1e-6;
    
    fmodel_d = fmodel;
    cmodel_d = cmodel;

    x = x0;
    pred = 0;
    [~, d] = l1_criticality_measure_and_descent_direction(fmodel_d, ...
                                                      cmodel_d, x, ...
                                                      mu, epsilon, ...
                                                      lb, ub);

    if ~isempty(find(x > ub | x < lb, 1))
        warning('cmg:out_of_bounds', 'Point already out of bounds');
    end

    for iter = 1:max_iters
        
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
        x_prev = x;
        pred_prev = pred;
        x = x + t*d;
        x(lower_bound_hits) = lb(lower_bound_hits);
        x(upper_bound_hits) = ub(upper_bound_hits);
        
        s = x - x0;
        pred = predict_descent(fmodel, cmodel, x - x0, mu);
        if pred < pred_prev
            % Rollback
            x = x_prev;
            break
        end
        % Break if step long enough
        if radius - norm(s, inf) < 0.05*radius
            break
        end
        fmodel_d = shift_model(fmodel, s);
        for k = 1:length(cmodel)
            cmodel_d(k) = shift_model(cmodel(k), s);
        end
        try
        [sigma, pseudo_steepest_descent] = l1_criticality_measure_and_descent_direction(fmodel_d, ...
                                                          cmodel_d, ...
                                                          x, mu, ...
                                                          epsilon, ...
                                                          lb, ub);
        catch teste
            rethrow(teste);
        end
        
        if ~isempty(find(upper_bound_hits | lower_bound_hits, 1))
            d = pseudo_steepest_descent; % Restarting
        else
            % Compute average Hessian
            bpoint_prev = 0;
            H = zeros(dim);
            for k = 1:length(brpoints_crossed)
                bpoint = brpoints_crossed(k);
                interval_length = (bpoint - bpoint_prev);
                mid_interval = 0.5*interval_length;
                relative_mid_point = (x - x0) + (bpoint + mid_interval)*d;
                
                H = H + (interval_length/t)*l1_hessian(fmodel, cmodel, ...
                                                       mu, relative_mid_point);
                bpoint_prev = bpoint;
            end
            
            % Conjugate direction
            dHd = d'*H*d;
            if norm(dHd) < tol_h
                break
            end
            beta = -((pseudo_steepest_descent'*H*d)/dHd);
            %beta = 0
            d = pseudo_steepest_descent + beta*d;
            
        end
    end
end
