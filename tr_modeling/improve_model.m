function model = improve_model(model, ff, options)
% IMPROVE_MODEL improves model to be lambda-poised (for some
% lambda)


% pivot_threshold is constant through the algorithm. Ensures we can
% obtain lambda-poised models when needed
tol_pivot = options.pivot_threshold;

n_interpolating_functions = length(ff);

while true
    % Recomplete interpolation set and calculate new model
    [model, exitflag] = complete_interpolation_set(model, ff, options);
    if exitflag >= 0 || model.radius < options.tol_radius
        break
    else
        model.radius = 0.5*model.radius;
    end
end

radius = model.radius;

smallest_pivot =  min(model.pivot_absvalues);
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

    %%
    while true
        % Perform the actual completion of interpolation set
        [points_scaled, indices, pivot_absvalues, pivot_polynomials] = ...
                change_interpolation_set_lu(points_scaled, basis, tol_pivot, true);
        choose_pivot = find(pivot_absvalues < tol_pivot, 1);
        if ~isempty(choose_pivot)
            chosen_poly = pivot_polynomials(choose_pivot);
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            smaller_radius = min(0.75, radius/scale_factor_x);
            for attempt = 1:4
                if attempt == 3
                    % Minimize inside TR
                    [new_point, new_value] = ...
                        minimize_polynomial(chosen_poly, smaller_radius);
                    if abs(new_value) < tol_pivot
                        % Not worth evaluating f
                        continue
                    end
                elseif attempt == 4
                    % Maximize inside TR
                    [new_point, new_value] = ...
                        minimize_polynomial(multiply_p(chosen_poly, ...
                                                       -1), smaller_radius);
                    if abs(new_value) < tol_pivot
                        % Not worth evaluating f
                        continue
                    end
                elseif attempt == 1
                    % Minimize in larger region
                    [new_point, new_value] = minimize_polynomial(chosen_poly, 1);
                    if abs(new_value) < tol_pivot
                        % Not worth evaluating f
                        continue
                    end
                else
                    % Maximize in larger region
                    [new_point, new_value] = ...
                        minimize_polynomial(multiply_p(chosen_poly, -1), 1);
                    if abs(new_value) < tol_pivot
                        error('cmg:pivot_not_found', 'Could not find pivot');
                    end
                end
                % Good pivot value, worth evaluating df function
                new_point_abs = scale_factor_x*new_point + center;
                new_fvalues = zeros(n_interpolating_functions, 1);
                choose_other_point = false;
                for nf = 1:n_interpolating_functions
                    new_fvalues(nf, 1) = ff{nf}(new_point_abs);
                    if isinf(new_fvalues(nf, 1)) || isnan(new_fvalues(nf, 1))
                        if attempt < 4
                            choose_other_point = true;
                            break
                        else
                            % Error
                            error('cmg:bad_fvalue', 'Bad f value');
                            return
                            % Or maybe I should just fail to add the point
                            % or just choose another polynomial
                        end
                    end
                end
                if choose_other_point
                    continue
                else
                    points_scaled(:, end+1) = new_point;
                    points_abs(:, end+1) = new_point_abs;
                    fvalues(:, end+1) = new_fvalues;
                    break
                end
            end        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            break
        end
    end
    smallest_pivot = min(pivot_absvalues);
    %%
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
    model.pivot_absvalues = pivot_absvalues;
    model.pivot_polynomials = pivot_polynomials;

end

%%

model.poised_radius = radius;
model.poised_center = model.points(:, 1);
model.smallest_pivot = smallest_pivot;



end

