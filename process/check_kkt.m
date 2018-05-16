function [status, lgrad] = check_kkt(f, c, x, tol_c, tol_g)

    dim = size(x, 1);
    nconstraints = length(c);
    
    [fval, fg] = f(x);
    G = zeros(dim, nconstraints);
    cvals = zeros(nconstraints, 1);
    for k = 1:nconstraints
        [cvals(k), G(:, k)] = c{k}(x);
    end
    % Start with 'more active':
    [~, ind] = sort(abs(cvals));
    cvals = cvals(ind);
    G = G(:, ind);

    active_constraints = false(1, nconstraints);
    lgrad = inf;
    status = false;
    for k = 1:nconstraints
        
        if cvals(k) < tol_c
            active_constraints(k) = true;
        else
            break
        end
        
        multipliers = -(G(:, active_constraints)\fg);

        lgrad = G(:, active_constraints)*multipliers + fg;

        if ~sum(multipliers < 0) && norm(lgrad) < tol_g
            status = true;
            break
        end
    end

end