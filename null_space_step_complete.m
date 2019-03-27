function [s, status] = null_space_step_complete(fmodel, cmodel, mu, ...
                                                x0, ind_qr, Q, R, ...
                                                radius, lb, ub, multipliers)

    tol_l = 1e-4;                    
    [dim, r_cols] = size(R);
    N = Q(:, r_cols+1:end);
    kappa_easy = 0.1;
    
    pgrad_center = l1_pseudo_gradient_new(fmodel, cmodel, mu);
    pH = l1_pseudo_hessian(fmodel, cmodel, mu, ind_qr, multipliers);

    finish = false;
    
    % Safe default
    s = zeros(dim, 1);
    
    bl_active = false(dim, 1);
    bu_active = false(dim, 1);
    for tries = 1:2*dim
        r_cols = size(R, 2);
        N = Q(:, r_cols+1:end);

        pgrad = pgrad_center + pH*s;
        g_red = N'*pgrad;
        H_red = N'*pH*N;

%          try
%              % This should be done in the step calculation
%              chol(H_red);
%              status = true;
%          catch
%              Dv = eig(H_red);
%              va = min(Dv);
%              if va < -tol_l
%                  % Will avoid calculating the step
%                  status = false;
%                  break
%              else
%                  status = true;
%              end
%          end
        radius_icb = radius;
        while true
            try
                [d_red, d_change] = ms_step(H_red, g_red, radius_icb);
            catch myerror
                rethrow(myerror);
            end
            d = N*d_red;
            d(bl_active | bu_active) = ...
                zeros(sum(bl_active | bu_active), 1);
            x = project_to_bounds(x0 + s, lb, ub);
            [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, d);
            if tmax_bounds >= 1 ...
                    && norm(s + d) < radius*(1 + kappa_easy)
                % Successfully computed step
                s = s + d;
                finish = true;
                break
            elseif tmax_bounds < 1
                bounded_change = tmax_bounds*(g_red'*d_red) ...
                    + 0.5*tmax_bounds^2*(d_red'*H_red*d_red);
                if bounded_change < 0.5*d_change ...
                        || (bounded_change <= 0 ...
                            && radius_icb < kappa_easy*radius)
                    % Decrease OK. Use this direction
                    bl_active_new = t_lb == tmax_bounds;
                    bu_active_new = t_ub == tmax_bounds;
                    [Q, R, bl_included, bu_included] = ...
                        include_bounds_gradients(Q, R, bl_active_new, ...
                                                 bu_active_new);

                    bl_active = bl_active | bl_included;
                    bu_active = bu_active | bu_included;

                    s = s + tmax_bounds*d;
                    break
                else
                    % Decrease not OK, recompute direction
                    radius_icb = 0.5*radius_icb;
                end
            else
                radius_icb = 0.5*radius_icb;
            end
        end
        if finish || norm(s) - radius >= kappa_easy
            break
        end
    end
end
