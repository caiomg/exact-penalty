function [x, pred] = try_to_make_activities_exact(fmodel, cmodel, mu, ...
                                                  ind_eactive, ...
                                                  tr_center, xh, ...
                                                  radius, lb, ub)
% TRY_TO_RECOVER_FEASIBILITY - 
%   

    xv = l1_feasibility_correction(cmodel, ind_eactive, tr_center, ...
                                   xh, radius, lb, ub);

    % Just compare descent obtained. Could include backtracking and all
    pred_h = predict_descent(fmodel, cmodel, xh - tr_center, mu);
    pred_v = predict_descent(fmodel, cmodel, xv - tr_center, mu);
    
    if pred_v >= pred_h
        pred = pred_v;
        x = xv;
    else
        pred = pred_h;
        x = xh;
    end

end
