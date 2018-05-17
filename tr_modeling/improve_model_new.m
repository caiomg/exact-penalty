function model = improve_model_new(model, functions, options)

    
    % Blind attempt to complete interpolation set
    while true
        [model, point_added] = improve_model_incremental(model, ...
                                                         functions, ...
                                                         options);
        if ~point_added
            model.radius = 0.5*model.radius;
            if model.radius < options.tol_radius
                break
            end
        else
            if is_lambda_poised(model, options)
                break
            end
        end
        
    end

end