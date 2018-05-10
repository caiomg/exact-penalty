function model = update_model_new(model, functions, options)



n_interpolating_functions = length(functions);

xi_imp = options.pivot_imp;
tol_pivot = options.pivot_threshold;
% First discart points which are more than a radius away
radius = model.radius;
poised_radius_factor = options.poised_radius_factor;
scale_factor_x = model.scale_factor_x;
if isempty(scale_factor_x) || scale_factor_x > radius*poised_radius_factor
    scale_factor_x = radius*poised_radius_factor;
end
if scale_factor_x < radius
   scale_factor_x = radius*poised_radius_factor; 
end
model = discard_points(model, scale_factor_x);


new_values_added = 0;

points_abs = model.points;
center = points_abs(:, 1);
fvalues = model.fvalues;

[dimension, p_ini] = size(points_abs);

n_interpolating_functions = length(functions);

% Scaling points to be in a ball of radius 1 around origin
points_scaled = points_abs;
max_norm = 0;
for column = 1:p_ini
    points_scaled(:, column) = points_abs(:, column) - center;
    points_scaled(:, column) = points_scaled(:, column)/scale_factor_x;
end

%%
basis = natural_basis(dimension);
basis_size = length(basis);


[points_scaled, indices, smallest_pivot, pivot_polynomials, is_pivot] = ...
                    check_geometry_lu_preserving(points_scaled, basis, tol_pivot);

% Reordering information AND deleting unused points
points_scaled = points_scaled(:, indices);
fvalues = fvalues(:, indices);
points_abs = points_abs(:, indices);



%%
n_points = size(points_scaled, 2);
basis_quadratic_size = dimension*(dimension + 1)/2;
% Calculate coefficients of polynomial model:
M = zeros(n_points, basis_size);
for m = 1:n_points
    for n = 1:basis_size
        M(m, n) = evaluate_polynomial(basis(n), points_scaled(:, m));
    end
end


if sum(is_pivot(1:dimension+1)) < dimension + 1
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
interp_tol = 5e-6;
% if abs(approaches_ndiff) > approaches_tol
%    error();
% end
% if norm(M*coefficients_basis - fvalues') > interp_tol
%     error();
% end
% if norm(M*coefficients2 - fvalues') > interp_tol
%     error();
% end


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
model.pivots = smallest_pivot;




end