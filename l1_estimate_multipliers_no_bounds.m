function [multipliers, tol] = l1_estimate_multipliers_no_bounds(p_grad, ...
                                                      Q, R)
% L1_ESTIMATE_MULTIPLIERS_NO_BOUNDS - 
%   
    ut_option.UT = true;
    r_cols = size(R, 2);
    
    % Actual calculation
    multipliers = -linsolve(R(1:r_cols, 1:r_cols), ...
                            Q(:, 1:r_cols)'*p_grad, ut_option);
    remainder = -(Q*(R*multipliers) + p_grad);

    % Correction
    [correction, rc] = linsolve(R(1:r_cols, 1:r_cols), ...
                                Q(:, 1:r_cols)'*remainder, ut_option);

    % Final calculation
    multipliers = multipliers + correction;

    % Tolerance
    tol = 10*max(norm(correction), eps);

end
