function model = try_to_add_interpolation_point(model, new_point, ...
                                                new_point_fvals, ...
                                                functions, options)

    % Add new point and its value to list of known points
    % Important not to add at the begining (which is the center)
    n_functions = length(functions);
    n_values = length(new_point_fvals);
    if n_functions ~= n_values
        error();
    end
    n_points = size(model.points, 2);
    
    model.points(:, n_points+1) = new_point;
    for nf = 1:n_functions
        model.fvalues(nf, n_points+1) = new_point_fvals{nf};
    end

    % Model improvement algorithm will improve poisedness of model
    % new point may be added or not
    model = improve_model(model, functions, options);
    
end