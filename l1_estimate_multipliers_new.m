function [multipliers, lb_multipliers, ub_multipliers, tol] = ...
        l1_estimate_multipliers_new(p_grad, Q, R, lb_active, ub_active)
% L1_ESTIMATE_MULTIPLIERS_NEW - 
%   
    dim = size(p_grad, 1);
    if isempty(find(lb_active | ub_active, 1))
        [multipliers, tol] = l1_estimate_multipliers_no_bounds(p_grad, Q, R);
        lb_multipliers = zeros(dim, 1);
        ub_multipliers = zeros(dim, 1);
    else
       [multipliers, lb_multipliers, ub_multipliers, tol] = ...
           l1_estimate_multipliers_with_bounds(p_grad, Q, R, lb_active, ...
                                               ub_active);
    end

end
