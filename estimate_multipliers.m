function multipliers = estimate_multipliers(fmodel, cmodel, x, mu, ...
                                            is_eactive, lb, ub)
    % ESTIMATE_MULTIPLIERS - 
    %   

    dim = size(x, 1);
    if sum(is_eactive) == 0
        multipliers = [];
    else
        s0 = zeros(dim, 1);
        pg = l1_pseudo_gradient_general(fmodel, cmodel, mu, s0, is_eactive);

        A = [cmodel(is_eactive).g];

        %% Indirect computation
        tol_bounds = min(1e-5, 1e-3*(ub - lb));
        lb_nearly_active = x - lb < tol_bounds;
        ub_nearly_active = ub - x < tol_bounds;

        n_nearly_active_nl_constraints = sum(is_eactive);
        n_nearly_active_bounds = sum(lb_nearly_active) + sum(ub_nearly_active);

        I = eye(dim);
        A_ext = [A, -I(:, lb_nearly_active), I(:, ub_nearly_active)];



        H = A_ext'*A_ext;
        g = A_ext'*pg;

        l_lb = [-inf(n_nearly_active_nl_constraints, 1);
                zeros(n_nearly_active_bounds, 1)];
        lambda0 = zeros(size(l_lb));

        lambda = solve_quadratic_problem(H, g, 0, [], [], [], [], l_lb, [], lambda0);

        multipliers = lambda(1:n_nearly_active_nl_constraints);
        % I should take note of which constraints these are

    end
end
