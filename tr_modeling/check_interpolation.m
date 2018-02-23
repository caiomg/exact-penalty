function max_diff = check_interpolation(model)

% Remove center from all points
h = model.points;
n_points = size(model.points, 2);
for m = n_points:-1:1
    h(:, m) = h(:, m) - model.points(:, 1);
end

max_diff = -1;
n_functions = size(model.fvalues, 1);
for k = 1:n_functions
    [c, g, H] = get_model_matrices(model, k-1);
    for m = 1:n_points
        this_value = c + g'*h(:, m) + 0.5*(h(:, m)'*H*h(:, m));
        difference = abs(model.fvalues(k, m) - this_value);
        if difference > max_diff
            max_diff = difference;
        end
        if abs(difference) > max(10*eps(max(abs(model.fvalues(k, :)))), sqrt(eps))
            error('cmg:tr_interpolation_error', 'Interpolation error');
        end
    end
end
