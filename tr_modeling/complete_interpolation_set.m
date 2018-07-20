function [model, exitflag] = complete_interpolation_set(model, ff, bl, bu, options)
% COMPLETE_INTERPOLATION_SET -- ensures the model has points so
% that interpolation is possible

dim = size(model.points, 1);
if isempty(bl)
    bl = -inf(dim, 1);
end
if isempty(bu)
    bu = inf(dim, 1);
end

n_interpolating_functions = length(ff);

if strcmp(options.basis, 'dummy')
    model.scale_factor_x = 1;
    center = model.points(:, 1);
    [c, g, H] = ff{1}(center);
    model.model_polynomial = matrices_to_polynomial(c, g, H);

    for nf = n_interpolating_functions:-1:2
        [c, g, H] = ff{nf}(center);
        model.other_polynomials{nf-1} = matrices_to_polynomial(c, g, H);
    end
    exitflag = 0;
    return
end

    
tol_pivot = options.pivot_threshold/30;
radius_factor = options.poised_radius_factor;

radius = model.radius;
points_abs = model.points;
fvalues = model.fvalues;
keep_points = discard_points_only(points_abs, radius*radius_factor);

points_abs = points_abs(:, keep_points);
fvalues = fvalues(:, keep_points);
center = points_abs(:, 1);
basis = model.basis;

basis_size = length(basis);

[dimension, p_ini] = size(points_abs);


% Scaling points to be in a ball of radius 1 around origin
points_scaled = points_abs;
max_norm = 0;
for column = 1:p_ini
    points_scaled(:, column) = points_abs(:, column) - center;
    norm_current = norm(points_scaled(:, column), 2);
    max_norm = max(max_norm, norm_current);
end
if max_norm > 0
    scale_factor_x = max(radius, max_norm);
else
    scale_factor_x = radius;
end
for column = 1:p_ini
    points_scaled(:, column) = points_scaled(:, column)/scale_factor_x;
end
bl_scaled = (bl - center)/scale_factor_x;
bl_scaled = min(bl_scaled, 0);
bu_scaled = (bu - center)/scale_factor_x;
bu_scaled = max(bu_scaled, 0);
unshift_point = @(x) x*scale_factor_x + center;
tol_shift = 10*eps(max(1, max(abs(center))))/min(1, scale_factor_x);
if tol_shift > 1e-3
    tol_shift = 1e-3; % DANGER, still
end
larger_radius = 1 - tol_shift;

pivot_polynomials = basis;
start_i = 1;
       while true
            % Perform the actual completion of interpolation set
            [points_scaled, indices, pivot_absvalues, pivot_polynomials] = ...
                    change_interpolation_set_lu(points_scaled, pivot_polynomials, tol_pivot, true, start_i);
            % Reordering information
            points_abs = points_abs(:, indices);
            points_scaled = points_scaled(:, indices);
            fvalues = fvalues(:, indices);
            % Finding point to change
            pivot_index = find(pivot_absvalues < tol_pivot, 1);
            smaller_radius = min(1-2*tol_pivot, radius/scale_factor_x);
            if ~isempty(pivot_index)
                chosen_poly = pivot_polynomials(pivot_index);
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                max_attempts = 3;
                for attempt = 1:max_attempts
                    pivot_found = false;
                    if attempt == 1
                        % Minimize inside TR
                        [new_point_a, new_value_a] = ...
                            minimize_polynomial_with_bounds(chosen_poly, bl_scaled, bu_scaled, smaller_radius);
                        % Maximize inside TR
                        [new_point_b, new_value_b] = ...
                            minimize_polynomial_with_bounds(multiply_p(chosen_poly, ...
                                                           -1), bl_scaled, bu_scaled, smaller_radius);
                        if abs(new_value_a) >= abs(new_value_b)
                            new_point = new_point_a;
                            new_value = new_value_a;
                        else
                            new_point = new_point_b;
                            new_value = new_value_b;
                        end
                    elseif attempt == 2
                        % Minimize inside TR
                        [new_point_a, new_value_a] = ...
                            minimize_polynomial_with_bounds(chosen_poly, bl_scaled, bu_scaled, larger_radius);
                        % Maximize inside TR
                        [new_point_b, new_value_b] = ...
                            minimize_polynomial_with_bounds(multiply_p(chosen_poly, ...
                                                           -1), bl_scaled, bu_scaled, larger_radius);
                        if abs(new_value_a) >= abs(new_value_b)
                            new_point = new_point_a;
                            new_value = new_value_a;
                        else
                            new_point = new_point_b;
                            new_value = new_value_b;
                        end
                    elseif attempt == 3
                        [new_point, new_value] = find_other_point_with_bounds(chosen_poly, bl_scaled, bu_scaled, larger_radius);
                    end
                    if abs(new_value) >= tol_pivot
                        pivot_found = true;
                        % Good pivot value, worth evaluating df function
                        new_point_abs = unshift_point(new_point);
                        new_fvalues = zeros(n_interpolating_functions, 1);
                        for nf = 1:n_interpolating_functions
                            new_fvalues(nf, 1) = ff{nf}(new_point_abs);
                            if isinf(new_fvalues(nf, 1)) || isnan(new_fvalues(nf, 1))
                                pivot_found = false;
                                break
                            end
                        end
                    end
                    if pivot_found
                        points_scaled = [points_scaled(:, 1:pivot_index-1), new_point, points_scaled(:, pivot_index:end)];
                        points_abs = [points_abs(:, 1:pivot_index-1), new_point_abs, points_abs(:, pivot_index:end)];
                        fvalues = [fvalues(:, 1:pivot_index-1), new_fvalues, fvalues(:, pivot_index:end)];
                        pivot_absvalues(pivot_index) = abs(new_value);
                        start_i = pivot_index;
                        break
                    end
                end
                if ~pivot_found
                    warning('cmg:bad_fvalue', 'Bad f value');
                    break
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                break
            end
        end
smallest_pivot = min(pivot_absvalues);
if smallest_pivot > tol_pivot
    exitflag = 0;
else
    exitflag = -1;
    return
end

% Discarting points unused for interpolation
points_abs = points_abs(:, 1:basis_size);
points_scaled = points_scaled(:, 1:basis_size);
fvalues = fvalues(:, 1:basis_size);

% I should scale points again



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
lastwarn('');
uf = linsolve(L, P*fvalues', l_opts);
coefficients_basis = linsolve(U, uf, u_opts);
[~, warnid] = lastwarn();
if ~isempty(warnid)
    1;
end
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
model.pivot_absvalues = pivot_absvalues;
model.pivot_polynomials = pivot_polynomials;



end


