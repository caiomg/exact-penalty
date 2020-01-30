function [x_trial, pred, Lambda] = l1_trust_region_step(fmodel, cmodel, x, ...
                                                epsilon, Lambda, mu, ...
                                                radius, lb, ub)
% L1_TRUST_REGION_STEP - 
%   



    [sigma, d, is_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel, cmodel, ...
                                                     x, mu, epsilon, lb, ub);
    xh = descent_direction_one_pass(fmodel, cmodel, mu, x, d, ...
                                    radius, lb, ub);
    % x_cg = conjugate_gradient_new_measure(fmodel, cmodel, x, mu, ...
    %                                       epsilon, radius, lb, ub, d);
    [x_trial, pred] = try_to_make_activities_exact(fmodel, cmodel, mu, ...
                                                   is_eactive, ...
                                                   x, xh, radius, lb, ub);
    if sigma < Lambda
        
        % Just reduce constraint values
        [x_other] = try_to_make_activities_exact(fmodel, cmodel, mu, ...
                                                 is_eactive, x, x, ...
                                                 radius, lb, ub, false);
        s = x_other - x;
        fmodel_shifted = shift_model(fmodel, s);
        for k = 1:length(cmodel)
            cmodel_shifted(k) = shift_model(cmodel(k), s);
        end
        
        % I should re-identify constraints
        multipliers = estimate_multipliers(fmodel_shifted, cmodel_shifted, ...
                                           x, mu, is_eactive, lb, ub);
        
        % Test some conditions for other descent directions
        tol_multipliers = 10*eps(mu);
        if ~isempty(find(multipliers < 0 - tol_multipliers, 1))
            warning('cmg:multipliers_negative', ...
                    'Negative Lagrange multiplier');
        end
        if ~isempty(find(multipliers > mu + tol_multipliers, 1))
            warning('cmg:multipliers_high', ...
                    'Overly high Lagrange multiplier');
        end
        try
        xm = l1_step_with_multipliers(fmodel_shifted, cmodel_shifted, ...
                                      x_other, multipliers, is_eactive, ...
                                      mu, x, radius, lb, ub);
        catch aqui
            rethrow(aqui);
        end
        pred_xm = predict_descent(fmodel, cmodel, xm - x, mu);
        if pred_xm > pred
            pred = pred_xm;
            x_trial = xm;
        else
            Lambda = 0.5*Lambda;
        end
    end
    
end
