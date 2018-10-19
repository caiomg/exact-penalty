function [model, exitflag] = try_to_add_point(model, new_point, new_fvalues, funcs, bl, bu, options)

    point_added = false;
    
    if ~is_complete(model)
        % Add this point
        relative_pivot_threshold = options.add_threshold;
        [model, point_added] = add_point(model, new_point, new_fvalues, relative_pivot_threshold);
    end
    if ~point_added
        model.cached_points = [new_point, model.cached_points];
        model.cached_fvalues = [new_fvalues, model.cached_fvalues];
        [model, exitflag] = ensure_improvement(model, funcs, bl, bu, options);
    else
        exitflag = 1;
    end
end