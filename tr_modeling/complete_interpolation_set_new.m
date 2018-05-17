function model = complete_interpolation_set_new(model, functions, options)

    
    % Blind attempt to complete interpolation set
    while true
        [model, point_added] = improve_model_incremental(model, ...
                                                         functions, ...
                                                         options);
        remaining = length(model.basis) - size(model.points, 2);
        if remaining <= 0
            break
        end
        if ~point_added
            model.radius = 0.5*model.radius;
            if model.radius < options.tol_radius
                break
            end
        end
    end

end