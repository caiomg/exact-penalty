function result = is_lambda_poised(model, options)
% IS_LAMBDA_POISED tests wether a model is lambda-poised for the
% given options
%
% Important: it assumes there is a finite radius_max defined.

% pivot_threshold defines how well-poised we demmand a model to be
xi = options.pivot_threshold;
poised_radius_factor = options.poised_radius_factor;
[dimension, n_points] = size(model.points);
tol_r = 1e-5;

result = false;
% Enough points
if n_points >= length(model.basis)
    % Passes pivot threshold
    if model.smallest_pivot >= xi
        % Radius not smaller than what is considered for poisedness
        if model.radius - model.scale_factor_x/poised_radius_factor >= -tol_r
            % Same center
            if norm(model.points(:, 1) - model.poised_center) == 0
                result = true;
            end
        end
    end
end

end