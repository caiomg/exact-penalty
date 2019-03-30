function [multipliers, tol, bl_mult, bu_mult] = ...
        l1_estimate_multipliers(fmodel, cmodel, mu, ind_eactive, Q, ...
                                R, x, bl, bu)

    if isempty(bl)
        bl = -inf(size(x));
    end
    if isempty(bu)
        bu = inf(size(x));
    end
    ut_option.UT = true;
    dimension = size(fmodel.g, 1);
    
    tol_rem = 1e-4;
    tol_mult = 1e-4;
    tol_ort = 1e-5;
    tol_con = 1e-10;
    
    n_constraints = size(R, 2);
    
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, cmodel, ...
                                         ind_eactive, true);

    [Q, R, lower_included, upper_included] = ...
        detect_and_include_active_bounds(Q, R, x, -pseudo_gradient, bl, bu, tol_con);

    cols_r = size(R, 2);

    % Actual multiplier estimation
    lastwarn(''); % Reseting warnings
    % First estimate
    multipliers = -linsolve(R(1:cols_r, 1:cols_r), (Q(:, 1:cols_r)'*pseudo_gradient), ut_option);
    [~, warnid] = lastwarn();
    if strcmp('MATLAB:nearlySingularMatrix', warnid)
        1; % Breakpoint for debugging
    end
    remainder = -((Q*R)*multipliers + pseudo_gradient);

    % Correction
    correction = linsolve(R(1:cols_r, 1:cols_r), Q(:, 1:cols_r)'*remainder, ut_option);
    % Final calculation
    multipliers = multipliers + correction;
    % Tolerance
    tol = 10*max(norm(correction), eps);
    
    bl_mult = zeros(size(bl));
    bu_mult = zeros(size(bu));

    n_lower_included = sum(lower_included);
    if n_lower_included > 0
        bl_mult(lower_included) = multipliers(n_constraints+1: ...
                                              n_constraints+ ...
                                              n_lower_included);
    end
    n_upper_included = sum(upper_included);
    if n_upper_included > 0
        bu_mult(upper_included) = multipliers(n_constraints + ...
                                              n_lower_included + 1 ...
                                              :n_constraints + ...
                                              n_lower_included + ...
                                              n_upper_included);
    end
    multipliers = multipliers(1:n_constraints);
    

end