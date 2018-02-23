function [points, indices, smallest_pivot, pivot_polynomials] = ...
                    change_interpolation_set_lu(points, basis, tol_pivot)
% COMPLETE_INTERPOLATION_SET - completes set of interpolation points
% through LU factorization

if nargin < 3 || isempty(tol_pivot)
    tol_pivot = sqrt(eps);
end
polynomials = basis;
n_polynomials = length(polynomials);
[dimension, n_points] = size(points);
% indices variable keeps track  of the order in the set of points
indices = 1:n_points;
smallest_pivot = Inf;
for pol = 1:n_polynomials
    pol_current = polynomials(pol);
    % Iterate through known points looking for a good pivot
    max_abs = 0;
    max_abs_idx = -1;
    for k = pol:n_points
        value = evaluate_polynomial(pol_current, points(:, indices(k)));
        if abs(value) > abs(max_abs)
            max_abs = value;
            max_abs_idx = k;
        end
    end
    if abs(max_abs) > tol_pivot
        % Known point satisfies pivot tolerance -- switch points
        indices([pol, max_abs_idx]) = indices([max_abs_idx, pol]);
    else
        % No known point satisfies tolerances
        % New point is calculated and added to the end of list of points
        new_point = find_other_point(pol_current);
        points(:, end+1) = new_point;
        n_points = n_points + 1;
        indices(end+1) = n_points;
        % Switch points
        indices([pol, n_points]) = indices([n_points, pol]);
    end
    value = evaluate_polynomial(pol_current, points(:, indices(pol)));
    if abs(value) < smallest_pivot
        smallest_pivot = abs(value);
    end
    % Make all remaining polynomials zero on current pivot
    for k = pol+1:n_polynomials
        polynomials(k) = zero_at_point(polynomials(k), polynomials(pol), ...
                                       points(:, indices(pol)));
    end
end
pivot_polynomials = polynomials;

end


    
