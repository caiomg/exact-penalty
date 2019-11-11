function [x_trial, pred] = l1_trust_region_step(fmodel, cmodel, x, ...
                                                epsilon, Lambda, mu, ...
                                                radius, lb, ub)
% L1_TRUST_REGION_STEP - 
%   



    [sigma, d, ind_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel, cmodel, ...
                                                     x, mu, epsilon, lb, ub);
    if sigma > Lambda
        xh = descent_direction_one_pass(fmodel, cmodel, mu, x, d, ...
                                        radius, lb, ub);
        [x_trial, pred] = try_to_make_activities_exact(fmodel, cmodel, mu, ...
                                                       ind_eactive, ...
                                                       x, xh, radius, lb, ub);
    else
        multipliers = estimate_multipliers(fmodel, cmodel, x, mu, ...
                                           ind_eactive, lb, ub);
        
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
        
        if sigma > 0
        % In this test I will be ignoring the multipliers
        % for lack of time to build the step
        xh = descent_direction_one_pass(fmodel, cmodel, mu, x, d, ...
                                        radius, lb, ub);
        [x_trial, pred] = try_to_make_activities_exact(fmodel, cmodel, mu, ...
                                                       ind_eactive, ...
                                                       x, xh, radius, lb, ub);
        else
            x_trial = x;
            pred = 0;
        end
    end
end
