function [x, pred] = try_to_make_activities_exact(fmodel, cmodel, mu, ...
                                                  ind_eactive, ...
                                                  tr_center, xh, ...
                                                  radius, lb, ub, guard_descent)
% TRY_TO_RECOVER_FEASIBILITY - 
%   

    if nargin < 10
        guard_descent = true;
    end
        
    xv = l1_feasibility_correction(cmodel, ind_eactive, tr_center, ...
                                   xh, radius, lb, ub);

    pred_v = predict_descent(fmodel, cmodel, xv - tr_center, mu);
    
    % Just compare descent obtained. Could include backtracking and all
    if guard_descent
        pred_h = predict_descent(fmodel, cmodel, xh - tr_center, mu);
        if pred_v >= pred_h
            pred = pred_v;
            x = xv;
        else
            pred = pred_h;
            x = xh;
        end
    else
        x = xv;
        pred = pred_v;
    end

end
