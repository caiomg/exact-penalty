function d_proj = projected_direction(x, d, Q, R, lb, ub)

    if isempty(lb)
        lb = -inf(size(x));
    end
    if isempty(ub)
        ub = inf(size(x));
    end

    dim = size(x, 1);
    tol_con = 1e-6;

    [Q, R, bl_included, bu_included] = ...
        detect_and_include_active_bounds(Q, R, x, d, lb, ub, tol_con);

    cols_r = size(R, 2);
    N = Q(:, cols_r+1:end);

    if isempty(N)
        d_proj = zeros(dim, 1);
    else
        d_proj = N*(N'*d);
    end

end

