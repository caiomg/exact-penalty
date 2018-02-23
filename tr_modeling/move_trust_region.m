function model = move_trust_region(model, new_center, new_center_fvals, ...
                                   functions, options)
% MOVE_TRUST_REGION accepts a TR step. The new point is the
% center of the new trust-region.

    % Adding the point to the beggining of the list of points it will be
    % considered the center. Also, it is associated with the
    % monomial 1 (constant).
    model.points = [new_center, model.points];

    % Adding the corresponding known function value (not to be calculated
    % again)
    n_functions = length(functions);
    n_values = length(new_center_fvals);
    if n_functions ~= n_values
        error();
    end
    model.fvalues(:, 2:end+1) = model.fvalues;
    for nf = 1:n_functions
        model.fvalues(nf, 1) = new_center_fvals(nf);
    end
    
    % The model will be updated with this new point. In our pivotal
    % approach, other points may end substituted.
    model = complete_interpolation_set(model, functions, options);
end