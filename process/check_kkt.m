function [status, lgrad] = check_kkt(f, c, x, bl, bu, tol_c, tol_g)

    if isempty(bl)
        bl = -inf(size(x));
    end
    if isempty(bu)
        bu = inf(size(x));
    end

    dim = size(x, 1);
    nconstraints = length(c);
    
    [fval, fg] = f(x);
    G = zeros(dim, nconstraints);
    cvals = zeros(nconstraints, 1);
    for k = 1:nconstraints
        [cvals(k), G(:, k)] = c{k}(x);
    end
    lcvals = bl - x;
    ucvals = x - bu;
    I = eye(dim);
    cvals = [cvals; lcvals; ucvals];
    G = [G, -I, I];

    % Start with 'more active':
    [~, ind] = sort(abs(cvals));
    cvals = cvals(ind);
    G = G(:, ind);
    max_cvals = max(cvals);

    active_constraints = false(1, nconstraints);
    lgrad = fg;
    status = false;
    fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                   'Algorithm', 'active-set', ...
                                   'SpecifyObjectiveGradient', false);
    for k = 1:nconstraints+2*dim
        
        if abs(cvals(k)) < tol_c
            active_constraints(k) = true;
        else
            break
        end
        
        m0 = -(G(:, active_constraints)\fg);
        m0_z = isinf(m0) | isnan(m0);
        m0(m0_z) = zeros(sum(m0_z), 1);

        lgrad_diff = @(m) G(:, active_constraints)*m + fg;
        n_lgrad = @(m) 0.5*norm(lgrad_diff(m))^2;

        multipliers = fmincon(n_lgrad, m0, [], [], [], [], ...
                            zeros(k, 1), [], [], fmincon_options);

        lgrad = G(:, active_constraints)*multipliers + fg;
        if norm(lgrad) < tol_g && max_cvals < tol_c
            break
        end
    end
    if norm(lgrad) < tol_g && max_cvals < tol_c
        status = true;
    end

end