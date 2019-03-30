function [Q, R, bl_included, bu_included] = ...
    detect_and_include_active_bounds(Q, R, x, d, lb, ub, tol_con)

    dim = size(x, 1);
    bl_included = false(dim, 1);
    bu_included = false(dim, 1);
    for k = 1:2*dim
        active_bounds = bl_included | bu_included;
        cols_r = size(R, 2);
        N = Q(:, cols_r+1:end);
        N(active_bounds, :) = zeros(sum(active_bounds), dim - cols_r);
        if cols_r == dim
            % All bounds included
            break
        end
        bl_newly_included = false(dim, 1);
        bu_newly_included = false(dim, 1);
        % Test bounds
        d_proj = N*(N'*d);
        [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, d_proj);
        if tmax_bounds*norm(d_proj) < tol_con
            newly_active_lb = t_lb == tmax_bounds;
            newly_active_ub = t_ub == tmax_bounds;
            max_product = 0;
            for m = find(newly_active_lb)'
                this_product = norm(N(m, :));
                if this_product > max_product
                    max_product = this_product;
                    max_grad = zeros(dim, 1);
                    max_grad(m) = -1;
                end
            end
            for m = find(newly_active_ub)'
                this_product = norm(N(m, :));
                if this_product > max_product
                    max_product = this_product;
                    max_grad = zeros(dim, 1);
                    max_grad(m) = 1;
                end
            end
            if max_product > 0
                [Q, R] = qrinsert(Q, R, cols_r + 1, max_grad);
                [val, pos] = max(max_grad);
                if val > 0
                    bu_newly_included(pos) = true;
                else
                    [~, pos] = min(max_grad);
                    bl_newly_included(pos) = true;
                end
            end
        end
        if ~isempty(find(bl_newly_included | bu_newly_included, 1))
            bl_included = bl_included | bl_newly_included;
            bu_included = bu_included | bu_newly_included;
        else
            break
        end
    end

end