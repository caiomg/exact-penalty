function [lb_active, ub_active] = active_bounds(x, d, lb, ub, tol_x, ...
                                                tol_d)

    if nargin < 5 || isempty(tol_x)
        tol_x = 1e2*eps(x);
    end
    if nargin < 6 || isempty(tol_d)
        tol_d = 1e3*eps(1);
    end

    if isempty(lb)
        lb_active = false(size(x));
    else
        lb_dist = x - lb;
        tl = -lb_dist./d;
        lb_active = lb_dist < 0 ...
            | (0 <= tl & (lb_dist < tol_x | tl < tol_d));
    end
    if isempty(ub)
        ub_active = false(size(x));
    else
        ub_dist = ub - x;
        tu = ub_dist./d;
        ub_active = ub_dist < 0 ...
            | (0 <= tu & (ub_dist < tol_x | tu < tol_d));
    end
    
end


    