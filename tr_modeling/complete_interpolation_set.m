function [model, smallest_pivot] = complete_interpolation_set(model, functions, options)
% COMPLETE_INTERPOLATION_SET -- ensures the model has points so
% that interpolation is possible

points_abs = model.points;
center = points_abs(:, 1);
basis = model.basis;
fvalues = model.fvalues;
radius = model.radius;

basis_size = length(basis);

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

% Perform the actual completion of interpolation set
[points_scaled, indices, smallest_pivot] = ...
        change_interpolation_set_lu(points_scaled, basis);

% Calculate new points in absolute coordinates
for k = p_ini+1:size(points_scaled, 2)
    points_abs(:, k) = scale_factor_x*points_scaled(:, k) + center;
end

% Calculate new function values
for k = p_ini+1:size(points_scaled, 2)
    for nf = 1:n_interpolating_functions
        fvalues(nf, k) = functions{nf}(points_abs(:, k));
    end
end

% Reordering information
points_abs = points_abs(:, indices);
points_scaled = points_scaled(:, indices);
fvalues = fvalues(:, indices);

% Discarting points unused for interpolation
points_abs = points_abs(:, 1:basis_size);
points_scaled = points_scaled(:, 1:basis_size);
fvalues = fvalues(:, 1:basis_size);

% Calculate coefficients of polynomial model:
M = zeros(basis_size);
for m = 1:basis_size
    for n = 1:basis_size
        M(m, n) = evaluate_polynomial(basis(n), points_scaled(:, m));
    end
end
[L, U, P] = lu(M);
l_opts.LT = true;
u_opts.UT = true;
uf = linsolve(L, P*fvalues', l_opts);
coefficients_basis = linsolve(U, uf, u_opts);

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
model.smallest_pivot = smallest_pivot;
model.poised_center = center;
model.poised_radius = scale_factor_x;



end


