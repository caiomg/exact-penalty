function x = descent_direction_one_pass(fmodel, cmodel, mu, x0, d, radius, lb, ub)

    x = project_to_bounds(x0, lb, ub);
    [t_lb, t_ub, tmax_bounds] = bounds_breakpoints(x, lb, ub, d);
    tmax_tr = infinity_tr_radius_breakpoint(d, radius);
    
    tmax = min(tmax_bounds, tmax_tr);
    try
        [t, brpoints_crossed] = line_search_cg(fmodel, cmodel, mu, d, tmax);
    catch this_error
        rethrow(this_error)
    end
    
    lower_bound_hits = t_lb == t;
    upper_bound_hits = t_ub == t;
    x = x + t*d;
    x(lower_bound_hits) = lb(lower_bound_hits);
    x(upper_bound_hits) = ub(upper_bound_hits);
    

end
