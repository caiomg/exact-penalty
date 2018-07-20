function [multipliers, tol, bl_mult, bu_mult] = l1_estimate_multipliers(fmodel, cmodel, mu, ind_eactive, Q, R, x, bl, bu)

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
    
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, cmodel, ...
                                         ind_eactive, true);

	% Test bounds
    l_active = x - bl <= tol_con & pseudo_gradient > 0;
    u_active = x - bu >= -tol_con & pseudo_gradient < 0;
    n_lower_active = sum(l_active);
    n_upper_active = sum(u_active);
    I = eye(dimension);
    Al = -I(:, l_active);
    Au = I(:, u_active);
    bounds_included = 0;
    lower_included = 0;
    upper_included = 0;
    for k = 1:n_lower_active
        this_grad = Al(:, k);
        norm_n = norm((Q*R)*(R\(Q'*this_grad)) - this_grad, 1);
        if norm_n > tol_ort
            lower_included = lower_included + 1;
            bounds_included = bounds_included + 1;
            [Q, R] = qrinsert(Q, R, bounds_included, this_grad);
        else
            l_active(this_grad ~= 0) = false;
        end
    end
    for k = 1:n_upper_active
        this_grad = Au(:, k);
        norm_n = norm((Q*R)*(R\(Q'*this_grad)) - this_grad, 1);
        if norm_n > tol_ort
            upper_included = upper_included + 1;
            bounds_included = bounds_included + 1;
            [Q, R] = qrinsert(Q, R, bounds_included, this_grad);
        else
            u_active(this_grad ~= 0) = false;
        end
    end
    if lower_included ~= sum(l_active)
        error()
    end
    if upper_included ~= sum(u_active)
        error()
    end
    
    % Estimate multipliers
    cols_r = size(R, 2);
    
    lastwarn(''); % Reseting warnings
    % First estimate
    multipliers = -linsolve(R(1:cols_r, :), (Q(:, 1:cols_r)'*pseudo_gradient), ut_option);
    [~, warnid] = lastwarn();
    if strcmp('MATLAB:nearlySingularMatrix', warnid)
        1; % Breakpoint for debugging
    end
    remainder = -((Q*R)*multipliers + pseudo_gradient);

    % Correction
    correction = linsolve(R(1:cols_r, :), Q(:, 1:cols_r)'*remainder, ut_option);
    % Final calculation
    multipliers = multipliers + correction;
    % Tolerance
    tol = 10*max(norm(correction), eps);
    
    bl_mult = zeros(size(bl));
    bu_mult = zeros(size(bu));

    if lower_included > 0
        bl_mult(l_active) = multipliers(1:lower_included);
    end
    if upper_included > 0
        bu_mult(u_active) = multipliers(lower_included + 1:lower_included + upper_included);
    end
    multipliers = multipliers(lower_included + upper_included + 1:end);
    

end