function [multipliers, tol, bl_mult, bu_mult] = l1_estimate_multipliers(fmodel, cmodel, mu, ind_eactive, Q, R, N, x, bl, bu)

    if isempty(bl)
        bl = -inf(size(x0));
    end
    if isempty(bu)
        bu = inf(size(x0));
    end
    ut_option.UT = true;
    dimension = size(fmodel.g, 1);
    
    tol_rem = 1e-4;
    tol_mult = 1e-4;
    
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, cmodel, ...
                                         ind_eactive, true);

	% Test bounds
    l_active = x <= bl;
    u_active = x >= bu;
    n_lower_active = sum(l_active);
    n_upper_active = sum(u_active);
    I = eye(dimension);
    A1 = [I(:, l_active), -I(:, u_active)];
    A = [A1, Q*R];
    [Q, R] = qr(A);
    N = null(A');
    
    % Estimate multipliers
    rows_qr = size(R, 1) - size(N, 2);
    
    lastwarn(''); % Reseting warnings
    % First estimate
    multipliers = -linsolve(R(1:rows_qr, :), (Q(:, 1:rows_qr)'*pseudo_gradient), ut_option);
    [~, warnid] = lastwarn();
    if strcmp('MATLAB:nearlySingularMatrix', warnid)
        1; % Breakpoint for debugging
    end
    remainder = -((Q*R)*multipliers + pseudo_gradient);

    % Correction
    correction = linsolve(R(1:rows_qr, :), Q(:, 1:rows_qr)'*remainder, ut_option);
    % Final calculation
    multipliers = multipliers + correction;
    % Tolerance
    tol = 10*max(norm(correction), eps);
    
    bl_mult = zeros(size(bl));
    bu_mult = zeros(size(bu));

    if n_lower_active > 0
        bl_mult(l_active) = multipliers(1:n_lower_active);
    end
    if n_upper_active > 0
        bu_mult(u_active) = multipliers(n_lower_active + 1:n_lower_active + n_upper_active);
    end
    multipliers = multipliers(n_lower_active + n_upper_active + 1:end);
    

end