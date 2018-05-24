function [model, new_values_added] = ...
        complete_interpolation_set_incremental(model, functions, options)



n_interpolating_functions = length(functions);

tol_pivot = options.pivot_threshold;
% First discart points which are more than a radius away
radius = model.radius;
poised_radius_factor = 1;
scale_factor_x = model.scale_factor_x;
if isempty(scale_factor_x) || scale_factor_x > radius*poised_radius_factor
    scale_factor_x = radius*poised_radius_factor;
end
if scale_factor_x < radius
   scale_factor_x = radius*poised_radius_factor; 
end


points_abs = model.points;
center = points_abs(:, 1);
fvalues = model.fvalues;

[dimension, p_ini] = size(points_abs);
n_interpolating_functions = length(functions);

new_values_added = 0;

% Scaling points to be in a ball of radius 1 around origin
points_scaled = points_abs;
for column = 1:p_ini
    points_scaled(:, column) = points_abs(:, column) - center;
    points_scaled(:, column) = points_scaled(:, column)/scale_factor_x;
end

%%
basis = model.basis;
basis_size = length(basis);

[points_scaled, indices, smallest_pivot, pivot_polynomials, is_pivot] = ...
    check_geometry_lu(points_scaled, basis, tol_pivot);
not_pivot = ~is_pivot;

% Reordering information AND deleting unused points
newly_discarded = 1:p_ini;
newly_discarded(indices(indices <= p_ini)) = [];
% cached_points = [points_abs(:, newly_discarded), cached_points];
% cached_fvalues = [fvalues(:, newly_discarded), cached_fvalues];
points_scaled = points_scaled(:, indices);
fvalues = fvalues(:, indices);
points_abs = points_abs(:, indices);

old_npoints = size(points_abs, 2);
while sum(not_pivot)
    % At least one polynomial is not pivot
    [next_polynomial, index_p] = choose_pivot_polynomial(pivot_polynomials(not_pivot));
    indices = find(not_pivot);
    newindex = indices(index_p);
    tol_point = tol_pivot;

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
        new_fvalues = zeros(n_interpolating_functions, 1);
        choose_other_point = false;
        for nf = 1:n_interpolating_functions
            new_fvalues(nf, 1) = functions{nf}(new_point_abs);
            if isinf(new_fvalues(nf, 1)) || isnan(new_fvalues(nf, 1))
                if attempt < 4
                    choose_other_point = true;
                    break
                else
                    % Error
                    warning('cmg:bad_fvalue', 'Bad f value');
                    return
                    % Or maybe I should just fail to add the point
                    % or just choose another polynomial
                end
            end
        end
        if choose_other_point
            continue;
        end
        if newindex <= old_npoints
            fvalues(:, 1:old_npoints+1) = [fvalues(:, 1:newindex-1), new_fvalues, fvalues(:, newindex:end)];
        else
            fvalues(:, old_npoints+1) = new_fvalues;
        end
        new_values_added = new_values_added + 1;
        if newindex <= old_npoints
            points_scaled = [points_scaled(:, 1:newindex-1), new_point, points_scaled(:, newindex:end)];
            points_abs = [points_abs(:, 1:newindex-1), new_point_abs, points_abs(:, newindex:end)];
        else
            points_scaled(:, end+1) = new_point;
            points_abs(:, end+1) = new_point_abs;
        end
        if sum(not_pivot)
            % We have to add this point/polynomial as pivot
            is_pivot(newindex) = true;
        end
        break
    end
    [points_scaled, indices, smallest_pivot, pivot_polynomials, is_pivot] = ...
    check_geometry_lu(points_scaled, basis, tol_pivot);
    not_pivot = ~is_pivot;

    % Reordering information AND deleting unused points
    newly_discarded = 1:p_ini;
    newly_discarded(indices(indices <= p_ini)) = [];
    % cached_points = [points_abs(:, newly_discarded), cached_points];
    % cached_fvalues = [fvalues(:, newly_discarded), cached_fvalues];
    points_scaled = points_scaled(:, indices);
    fvalues = fvalues(:, indices);
    points_abs = points_abs(:, indices);    
    old_npoints = size(points_abs, 2);
end


n_points = size(points_scaled, 2);
% Calculate coefficients of polynomial model:
M = zeros(n_points, basis_size);
for m = 1:n_points
    for n = 1:basis_size
        M(m, n) = evaluate_polynomial(basis(n), points_scaled(:, m));
    end
end
if sum(is_pivot(1:dimension+1)) < dimension + 1 || sum(is_pivot(dimension+2:end)) < 1
    % Consider just the linear part
    r_rows = sum(is_pivot(1:dimension+1));
    Ml = M(1:r_rows, 1:dimension+1);
    [Q, R] = qr(Ml');
    linsolve_options.LT = true;
    coefficients_basis = zeros(basis_size, n_interpolating_functions);
    for nf = 1:n_interpolating_functions
        coefficients_basis(1:dimension+1, nf) = Q(:, 1:r_rows)*linsolve(R(1:r_rows, :)', fvalues(nf, 1:r_rows)', linsolve_options);
    end
else
    % Minimum norm hessian

    fmincon_options = optimoptions('fmincon', 'Display', 'off');
    linear_not_pivot_ind = 1:dimension + 1;
    linear_not_pivot_ind(is_pivot(1:dimension+1)) = [];
    quadratic_ind = dimension+2:basis_size;
    indices_to_minimize = [linear_not_pivot_ind, quadratic_ind];
    indices_to_minimize = [quadratic_ind];
    hnorm = @(x) 0.5*(x(indices_to_minimize)'*x(indices_to_minimize));
    for nf = 1:n_interpolating_functions
        % In the future I shoul apply a shift to remove the zero order term
        coefficients_basis(:, nf) = fmincon(hnorm, zeros(basis_size, 1), [], [], ...
                                            M, fvalues(nf, :)', ...
                                            [], [], [], fmincon_options);
    end
end
    



% % Second calculation
% [Q, R] = qr(M');
% linsolve_options.LT = true;
% % Least norm solution
% lns = Q(:, 1:n_points)*linsolve(R(1:n_points, :)', fvalues', linsolve_options);
% % we are assuming M has full row rank
% N = Q(:, n_points+1:basis_size);
% N2 = N(dimension+2:basis_size, :);
% %correction = (N2'*N2)\(-N2'*lns(dimension+2:end, :));
% correction = -(N2\lns(dimension+2:end, :));
% coefficients2 = lns + N*correction;
% 
% approaches_ndiff = hnorm(coefficients2) - hnorm(coefficients_basis);
% approaches_tol = 5e-6*max(1, hnorm(coefficients_basis));
% interp_tol = 5e-6;
% 
% % if abs(approaches_ndiff) > approaches_tol
% %    error();
% % end
% if norm(M*coefficients_basis - fvalues') > interp_tol
%     error();
% end
% % if norm(M*coefficients2 - fvalues') > interp_tol
% %     error();
% % end


% Calculate polynomial model
model_polynomial = polynomial_zero(dimension);
for m = 1:basis_size
   model_polynomial = add_p(model_polynomial, ...
                            multiply_p(basis(m), coefficients_basis(m, 1)));
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
model.poised_center = model.points(:, 1);
model.smallest_pivot = smallest_pivot;
model.poised_radius = radius;



end