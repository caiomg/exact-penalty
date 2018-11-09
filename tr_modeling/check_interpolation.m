function max_diff = check_interpolation(model)

% Tolerance
if model.radius < 1e-3 
    tol_1 = 100*eps;
else
    tol_1 = 10*eps;
end
tol_2 = 10*sqrt(eps);

% Remove shift center from all points
h = model.points_abs;
n_points = size(h, 2);
if n_points < model.tr_center
   1; 
end
for m = n_points:-1:1
    h(:, m) = h(:, m) - model.points_abs(:, model.tr_center);
end

n_functions = size(model.fvalues, 1);
max_diff = -ones(n_functions, 1);
tol = zeros(n_functions, 1);
for k = 1:n_functions
    [c, g, H] = get_model_matrices(model, k-1);
    tol(k) = max(tol_1*max(abs(model.fvalues(k, :))), tol_2);
    for m = 1:n_points
        this_value = c + g'*h(:, m) + 0.5*(h(:, m)'*H*h(:, m));
        difference = abs(model.fvalues(k, m) - this_value);
        if difference > max_diff(k)
            max_diff(k) = difference;
        end
    end
end
violated = (max_diff > tol);
max_diff = max(max_diff);
if ~isempty(find(violated, 1))
   warning('cmg:tr_interpolation_error', 'Interpolation error');
end

