function [measure, direction, removed] = l1_choose_search_direction_new(...
    Q, R, multipliers, pg, mu, x, lb, ub, tol_multipliers, threshold)

    if nargin < 10
        tol_multipliers = 0;
    end
    if nargin < 11
        threshold = inf;
    end
    
    n_constraints = length(multipliers);
    
    measure = 0;
    direction = [];
    Q_measure = Q;
    R_measure = R;
    removed = 0;

    for k = 1:n_constraints
        if measure >= threshold
            break
        end
        if multipliers(k) < -tol_multipliers ...
                || multipliers(k) - mu > tol_multipliers
            [Q_k, R_k] = qrdelete_fix(Q, R, k);
            con_grad = Q*R(:, k);
            if multipliers(k) > 0
                direction_k = con_grad;
                direction_k_proj = projected_direction(x, direction_k, Q_k, R_k, lb, ub);
                pgrad_k = pg + mu*con_grad;
            else
                direction_k = -con_grad;
                direction_k_proj = projected_direction(x, direction_k, Q_k, R_k, lb, ub);
                pgrad_k = pg;
            end
            pred_change = pgrad_k'*direction_k_proj;
            if pred_change < 0
                measure_k = sqrt(-pred_change);%/abs(multipliers(k));
            elseif pred_change == 0
                measure_k = 0;
            else
                measure_k = -1;
                warning('cmg:not_descent', ...
                        'Dropping step gave ascent %g', pred_change);
            end
            if measure_k > measure
                measure = measure_k;
                direction = direction_k_proj;
                Q_measure = Q_k;
                R_measure = R_k;
                removed = k;
            end
        end
    end
    Q = Q_measure;
    R = R_measure;
end


