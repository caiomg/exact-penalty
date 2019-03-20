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
    n_constraints = size(R, 2);
    lower_included = 0;
    upper_included = 0;
    for k = 1:n_lower_active
        this_grad = Al(:, k);
        cols_r = size(R, 2);
        norm_n = norm(Q(:, cols_r + 1:end)'*this_grad);
        if norm_n > tol_ort
            lower_included = lower_included + 1;
            [Q, R] = qrinsert(Q, R, cols_r+1, this_grad);
        else
            l_active(this_grad ~= 0) = false;
        end
    end
    for k = 1:n_upper_active
        this_grad = Au(:, k);
        cols_r = size(R, 2);
        norm_n = norm(Q(:, cols_r + 1:end)'*this_grad);
        if norm_n > tol_ort
            upper_included = upper_included + 1;
            [Q, R] = qrinsert(Q, R, cols_r+1, this_grad);
        else
            u_active(this_grad ~= 0) = false;
        end
    end
    cols_r = size(R, 2);
    assert(lower_included == sum(l_active));
    assert(upper_included == sum(u_active));
    assert(n_constraints + lower_included + upper_included == cols_r);

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

    if lower_included > 0
        bl_mult(l_active) = multipliers(n_constraints+1:n_constraints+lower_included);
    end
    if upper_included > 0
        bu_mult(u_active) = multipliers(n_constraints + lower_included + 1:n_constraints + lower_included + upper_included);
    end
    multipliers = multipliers(1:n_constraints);
    

end