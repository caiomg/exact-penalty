function [points, indices, pivot_absvalues, pivot_polynomials] = ...
                    change_interpolation_set_lu(points, basis, tol_pivot, ...
                    reorder, start_i, minimum_norm)
% COMPLETE_INTERPOLATION_SET - completes set of interpolation points
% through LU factorization

if nargin < 6
    minimum_norm = false;
end
if nargin < 5 || isempty(start_i)
    start_i = 1;
end

polynomials = basis;
n_polynomials = length(polynomials);
[dimension, n_points] = size(points);
pivot_absvalues = zeros(n_polynomials, 1);
% indices variable keeps track  of the order in the set of points
indices = 1:n_points;
point_pt = 1;
for pol = 1:start_i-1
    pol_current = polynomials(pol);
    value = evaluate_polynomial(pol_current, points(:, indices(point_pt)));
    pivot_absvalues(pol) = abs(value);
    point_pt = point_pt + 1;
end
for pol = start_i:min(n_polynomials, n_points)
    pol_current = polynomials(pol);
    % Iterate through known points looking for a good pivot
    max_abs = 0;
    max_abs_idx = -1;
    value = evaluate_polynomial(pol_current, points(:, indices(point_pt)));
    if abs(value) < tol_pivot || reorder
        for k = point_pt:n_points
            value = evaluate_polynomial(pol_current, points(:, indices(k)));
            if abs(value) > abs(max_abs)
                max_abs = abs(value);
                max_abs_idx = k;
            end
        end
        if abs(max_abs) > tol_pivot
            % Known point satisfies pivot tolerance -- switch points
            indices([point_pt, max_abs_idx]) = indices([max_abs_idx, point_pt]);
        else
            % No known point satisfies tolerances
            if minimum_norm
                pivot_absvalues(pol) = -0.5*tol_pivot;
                continue
            else
                break
            end
        end
    end
    value = evaluate_polynomial(pol_current, points(:, indices(pol)));
    pivot_absvalues(pol) = abs(value);
    % Make all remaining polynomials zero on current pivot
    for k = pol+1:n_polynomials
        polynomials(k) = zero_at_point(polynomials(k), polynomials(pol), ...
                                       points(:, indices(point_pt)));
    end
    point_pt = point_pt + 1;
end
pivot_polynomials = polynomials;

end


    
