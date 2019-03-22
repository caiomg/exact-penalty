function [Q, R, bl_included, bu_included] = ...
        include_bounds_gradients(Q, R, include_bl, include_bu)


    tol_indep = 1e-5;
    
    [dim, r_cols] = size(R);

    gradients_ub = eye(dim);
    gradients_lb = -gradients_ub;

    bl_included = false(dim, 1);
    bu_included = false(dim, 1);
    for k = 1:dim
        if include_bl(k)
            this_gradient = gradients_lb(:, k);
            independent = Q(:, r_cols+1:end)'*this_gradient;
            if norm(independent) > tol_indep
                [Q, R] = qrinsert(Q, R, r_cols + 1, this_gradient);
                r_cols = r_cols + 1;
                bl_included(k) = true;
            end
        end
        if include_bu(k)
            this_gradient = gradients_ub(:, k);
            independent = Q(:, r_cols+1:end)'*this_gradient;
            if norm(independent) > tol_indep
                [Q, R] = qrinsert(Q, R, r_cols + 1, this_gradient);
                r_cols = r_cols + 1;
                bu_included(k) = true;
            end
        end
    end

    
end