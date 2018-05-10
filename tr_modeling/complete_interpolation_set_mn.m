function model = complete_interpolation_set_mn(model, functions, options)

    
    [dimension, n_points] = size(model.points);
    % Blind attempt to complete interpolation set
    while ~is_lambda_poised(model, options)
        [model, point_added] = improve_model_incremental(model, ...
                                                         functions, ...
                                                         options);
        if ~point_added
            error('cmg:no_point_added', 'No point added');
        end
    end

end