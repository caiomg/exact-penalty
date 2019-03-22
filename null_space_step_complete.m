function [s, status] = null_space_step_complete(fmodel, cmodel, mu, ...
                                                x0, ind_qr, Q, R, ...
                                                radius, lb, ub, multipliers)

    tol_l = 1e-4;                    
    [dim, r_cols] = size(R);
    N = Q(:, r_cols+1:end);
    
    pgrad = l1_pseudo_gradient_new(fmodel, cmodel, mu);
    pH = l1_pseudo_hessian(fmodel, cmodel, mu, ind_qr, multipliers);
    H_red = N'*pH*N;


    
    % Safe default
    s = zeros(dim, 1);
    
    bl_active = false(dim, 1);
    bu_active = false(dim, 1);
    for tries = 1:2*dim
        r_cols = size(R, 2);
        N = Q(:, r_cols+1:end);
        
        g_red = N'*pgrad;
        H_red = N'*pH*N;

        try
            % This should be done in the step calculation
            chol(H_red);
            status = true;
        catch
            Dv = eig(H_red);
            va = min(Dv);
            if va < -tol_l
                % Will avoid calculating the step
                status = false;
                break
            else
                status = true;
            end
        end

        s_red = ms_step(H_red, g_red, radius);
        s = N*s_red;

        s(bl_active | bu_active) = ...
            zeros(sum(bl_active | bu_active), 1);
        
        % Identify activated bounds
        [bl_active_x, bu_active_x] = active_bounds(x0, s, lb, ub);
        bl_active_new = bl_active_x & ~bl_active;
        bu_active_new = bu_active_x & ~bu_active;
        if isempty(find(bl_active_new | bu_active_new, 1))
            break
        else
            [Q, R, bl_included, bu_included] = ...
                include_bounds_gradients(Q, R, bl_active_new, ...
                                         bu_active_new);
            bl_active = bl_active | bl_included;
            bu_active = bu_active | bu_included;
        end
    end
end
