function [model, new_values_added] = improve_model_incremental(model, functions, options)



n_interpolating_functions = length(functions);

% First discart points which are more than a radius away
radius = model.radius;
discard_tol = 1.05;
poised_radius_factor = options.poised_radius_factor*discard_tol;
model = discard_points(model, poised_radius_factor*radius);

xi_imp = options.pivot_imp;
tol_pivot = options.pivot_threshold;
new_values_added = 0;

points_abs = model.points;
center = points_abs(:, 1);
fvalues = model.fvalues;
radius = model.radius;

[dimension, p_ini] = size(points_abs);

n_interpolating_functions = length(functions);

% Scaling points to be in a ball of radius 1 around origin
points_scaled = points_abs;
max_norm = 0;
for column = 1:p_ini
    points_scaled(:, column) = points_abs(:, column) - center;
    norm_current = norm(points_scaled(:, column), 2);
    max_norm = max(max_norm, norm_current);
end
if max_norm > 0
    scale_factor_x = max_norm;
else
    scale_factor_x = radius;
end
for column = 1:p_ini
    points_scaled(:, column) = points_scaled(:, column)/scale_factor_x;
end

%%
basis = model.basis;
basis_size = length(basis);

[points_scaled, indices, smallest_pivot, pivot_polynomials, is_pivot] = ...
    check_geometry_lu(points_scaled, basis, tol_pivot);

% Reordering information AND deleting unused points
points_scaled = points_scaled(:, indices);
fvalues = fvalues(:, indices);
points_abs = points_abs(:, indices);

if sum(~is_pivot)
    % At least one polynomial is not pivot
    newindex = size(points_abs, 2) + 1;
    next_polynomial = choose_pivot_polynomial(pivot_polynomials(~is_pivot));
    tol_point = tol_pivot;
else
    % All polynomials are pivots
    newindex = size(points_abs, 2);
    next_polynomial = pivot_polynomials(end);
    current_value = evaluate_polynomial(next_polynomial, points_scaled(:, newindex));
    tol_point = max(tol_pivot, xi_imp*abs(current_value));
end
% Add to the set

% Trying to find a point strictly inside the current
% trust-region
scaled_radius = radius/scale_factor_x;
for attempt = 1:4
    if attempt == 1
        % Minimize inside TR
        [new_point, new_value] = minimize_polynomial(next_polynomial, ...
                                                     scaled_radius);
        if abs(new_value) < tol_point
            % Not worth evaluating f
            continue
        end
    elseif attempt == 2
        % Maximize inside TR
        [new_point, new_value] = ...
            minimize_polynomial(multiply_p(next_polynomial, -1), ...
                                scaled_radius);
        if abs(new_value) < tol_point
            % Not worth evaluating f
            continue
        end
    elseif attempt == 3
        % Minimize in larger region
        [new_point, new_value] = minimize_polynomial(next_polynomial, 1);
        if abs(new_value) < tol_point
            % Not worth evaluating f
            continue
        end
    else
        % Maximize in larger region
        [new_point, new_value] = ...
            minimize_polynomial(multiply_p(next_polynomial, -1), 1);
        if abs(new_value) < tol_point
            % Not worth evaluating f
            continue
        end
    end
    % Good pivot value, worth evaluating df function
    new_point_abs = scale_factor_x*new_point + center;
    attempt_failed = false;
    for nf = 1:n_interpolating_functions
        fvalues(nf, newindex) = functions{nf}(new_point_abs);
        if isinf(fvalues(nf, newindex)) || isnan(fvalues(nf, newindex))
            if attempt < 4
                attempt_failed = true;
                break
            else
                % Error
                error('cmg:bad_fvalue', 'Bad f value');
                % Or maybe I should just fail to add the point
                % or just choose another polynomial
            end
        end
    end
    if attempt_failed
        continue
    end
    points_scaled(:, newindex) = new_point;
    points_abs(:, newindex) = new_point_abs;
    new_values_added = new_values_added + 1;
    break
end
    


%%
n_points = size(points_scaled, 2);
% Calculate coefficients of polynomial model:
M = zeros(n_points, basis_size);
for m = 1:n_points
    for n = 1:basis_size
        M(m, n) = evaluate_polynomial(basis(n), points_scaled(:, m));
    end
end


fmincon_options = optimoptions('fmincon', 'Display', 'off');
hnorm = @(x) 0.5*(x(dimension+2:basis_size)'*x(dimension+2:basis_size));
for nf = 1:n_interpolating_functions
    % In the future I shoul apply a shift to remove the zero order term
    coefficients_basis(:, nf) = M\fvalues(nf, :)';
end
if sum(isnan(coefficients_basis))
    1;
end

% Calculate polynomial model
model_polynomial = polynomial_zero(dimension);
try
for m = 1:basis_size
   model_polynomial = add_p(model_polynomial, ...
                            multiply_p(basis(m), coefficients_basis(m, 1)));
end
catch erro
    rethrow(erro);
end
    
for nf = n_interpolating_functions:-1:2
    c_polynomial = polynomial_zero(dimension);
    for m = 1:basis_size
       c_polynomial = add_p(c_polynomial, ...
                                multiply_p(basis(m), coefficients_basis(m, nf)));
    end
    other_polynomials{nf-1} = c_polynomial;
end

model.model_polynomial = model_polynomial;
if n_interpolating_functions > 1
    model.other_polynomials = other_polynomials;
end
model.points = points_abs;
model.fvalues = fvalues;
model.scale_factor_x = scale_factor_x;
model.poised_radius = radius;
model.poised_center = model.points(:, 1);
model.smallest_pivot = smallest_pivot;




end