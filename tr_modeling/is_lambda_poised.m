function result = is_lambda_poised(model, options)
% IS_LAMBDA_POISED tests wether a model is lambda-poised for the
% given options
%
% Important: it assumes there is a finite radius_max defined.

% pivot_threshold defines how well-poised we demmand a model to be
xi = options.pivot_threshold;

if strcmp(options.basis, 'dummy')
    result = true;
    return
end
result = false;
% Passes pivot threshold
if model.smallest_pivot >= xi
    % Larger radius
    if model.radius >= model.poised_radius
        % Same center
        if norm(model.points(:, 1) - model.poised_center) == 0
            result = true;
        end
    end
end

end