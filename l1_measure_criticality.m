function [q1, q2, p_grad, d, Q, R, ind_qr, multipliers, tol_multipliers] = ...
        l1_measure_criticality(fmodel, cmodel, mu, Q, R, ind_qr, x, lb, ...
                                ub, threshold)
% L1_MEASURE_CRITICALITY - 
%   

    p_grad = l1_pseudo_gradient(fmodel.g, mu, cmodel, ind_qr, true);

    pg_proj = projected_direction(x, -p_grad, Q, R, lb, ub);
    q1 = norm(pg_proj);
    
    if q1 < threshold
        [multipliers, tol_multipliers, bl_mult, bu_mult, Qm, Rm] = ...
            l1_estimate_multipliers(fmodel, cmodel, mu, ind_qr, Q, ...
                                    R, x, lb, ub);
        if ~isempty(find(multipliers < -tol_multipliers ...
                | multipliers - mu > tol_multipliers, 1))
            1;
        end
        [q2, d, Q_drop, R_drop, ind_qr_drop] = l1_choose_search_direction(...
            Qm, Rm, ind_qr, multipliers, p_grad, mu, x, lb, ub, tol_multipliers, ...
            threshold);

        if q2 > q1
            Q = Q_drop;
            R = R_drop;
            ind_qr = ind_qr_drop;
        end
    else
        q2 = 0;
        d = [];
        multipliers = [];
    end
    
end