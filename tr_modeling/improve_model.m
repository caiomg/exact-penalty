function model = improve_model(model, functions, options)
% IMPROVE_MODEL improves model to be lambda-poised (for some
% lambda)


% pivot_threshold is constant through the algorithm. Ensures we can
% obtain lambda-poised models when needed
tol_pivot = options.pivot_threshold;

n_interpolating_functions = length(functions);

% First discart points which are more than a radius away
radius = model.radius;
model = discart_points(model, radius);

% Recomplete interpolation set and calculate new model
[model, smallest_pivot] = complete_interpolation_set_incremental(model, functions, options);

% Check if pivot_threshold is already satisfied
if smallest_pivot < tol_pivot

    points_abs = model.points;
    center = points_abs(:, 1);
    basis = model.basis;
    fvalues = model.fvalues;

    basis_size = length(basis);
    [dimension, p_ini] = size(points_abs);

    % Scaling
    points_scaled = points_abs;
    scale_factor_x = model.scale_factor_x;
    for column = 1:p_ini
        points_scaled(:, column) = (points_abs(:, column) ...
                                    - center)/scale_factor_x;
    end

    [points_scaled, indices, smallest_pivot] = ...
        change_interpolation_set_lu(points_scaled, basis, tol_pivot);

    % Calculate new points in absolute coordinates
    for k = p_ini+1:size(points_scaled, 2)
        points_abs(:, k) = scale_factor_x*points_scaled(:, k) + center;
    end
    
    % Calculate new function values
    for k = p_ini+1:size(points_scaled, 2)
        for nf = 1:n_interpolating_functions
            fvalues(nf, k) = functions{nf}(points_abs(:, k));
            if isnan(fvalues(nf, k)) || isinf(fvalues(nf, k))
                error('cmg:bad_fvalue', 'Bad function value');
            end
        end
    end

    % Reordering information
    points_abs = points_abs(:, indices);
    points_scaled = points_scaled(:, indices);
    fvalues = fvalues(:, indices);

    % Discarting points not used for interpolation
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

end

model.poised_radius = radius;
model.poised_center = model.points(:, 1);
model.smallest_pivot = smallest_pivot;



end

