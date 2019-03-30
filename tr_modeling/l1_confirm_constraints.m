function [Q, R, ind_qr, crit_measure, multipliers, tol_multipliers, dropped_gradient] = ...
        l1_confirm_constraints(fmodel, cmodel, ind_eactive, ind_qr, ...
                               Q, R, mu, x, bl, bu, Lambda)

   
    while true
        pgradient = ...
            l1_pseudo_gradient(fmodel.g, mu, cmodel, ind_qr, true);
        crit_measure = ...
            l1_criticality_measure(x, pgradient, Q, R, bl, bu, []);
        if crit_measure >= Lambda
            % No use calculating multipliers
            multipliers = [];
            tol_multipliers = 0;
            break
        end
        [multipliers, tol_multipliers] = ...
            l1_estimate_multipliers(fmodel, cmodel, ...
                                    mu, ind_qr, Q, R, x, bl, bu);
        if isempty(find(multipliers < -tol_multipliers ...
                         | mu < multipliers - tol_multipliers, 1))
            % Multipliers already on correct range
            break
        else
            % Need to drop a constraint
            [Q, R, ind_qr, ind_eactive, newly_dropped, multiplier_dropped] = ...
                l1_drop_constraint(cmodel, Q, R, ind_qr, ind_eactive, mu, ...
                                   multipliers, tol_multipliers);
%             cols_r = size(R, 2);
%             N = Q(:, cols_r+1:end);
%             dropped_gradient = cmodel(newly_dropped).g;
%             sigma = max(sign(multiplier_dropped));
%             expected_gradient = pgradient + sigma*mu*dropped_gradient;
%             s_direction = sigma*(N*(N'*dropped_gradient));
%             predicted_reduction = expected_gradient'*s_direction;
        end
    end
end
