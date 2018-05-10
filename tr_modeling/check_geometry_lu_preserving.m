function [points, pt_order, smallest_pivot, pivot_polynomials, poly_is_pivot] = ...
                    check_geometry_lu_preserving(points, basis, tol_pivot)
% CHECK_GEOMETRY_LU - Uses LU factorization to ensure points are linearly
% independent with regards to the polynomials. The interpolation matrix may
% not have full column rank.
% May be used for discarding points or finding the best pivot polynomial to
% improve.

if nargin < 3 || isempty(tol_pivot)
    tol_pivot = sqrt(eps);
end
polynomials = basis;
n_polynomials = length(polynomials);
[dimension, n_points] = size(points);
% indices variable keeps track  of the order in the set of points
pt_order= 1:n_points;
all_polynomials = 1:n_polynomials;
poly_is_pivot = false(1, n_polynomials);
smallest_pivot = Inf;
point_pt = 1;
for pol = 1:n_polynomials
    pol_current = polynomials(pol);
    % Iterate through known points looking for a good pivot
    max_abs = 0;
    max_abs_idx = -1;
    value = evaluate_polynomial(pol_current, points(:, pt_order(point_pt)));
    if abs(value) > tol_pivot
        max_abs = value;
        max_abs_idx = point_pt;
    else
        min_dist = inf;
        for k = point_pt+1:n_points
            value = evaluate_polynomial(pol_current, points(:, pt_order(k)));
            if abs(value) > tol_pivot && min_dist > norm(points(:, pt_order(k)))
                max_abs = value;
                max_abs_idx = k;
                min_dist = norm(points(:, pt_order(k)));
            end
        end
    end
    if abs(max_abs) > tol_pivot
        % Known point satisfies pivot tolerance -- switch points
        pt_order([point_pt, max_abs_idx]) = pt_order([max_abs_idx, point_pt]);
        value = evaluate_polynomial(pol_current, points(:, pt_order(point_pt)));
        if abs(value) < smallest_pivot
            smallest_pivot = abs(value);
        end
        poly_is_pivot(pol) = true;
        %remaining_polynomials = all_polynomials(~poly_is_pivot);
        % Make all remaining polynomials zero on current pivot
        for k = pol+1:n_polynomials
            polynomials(k) = zero_at_point(polynomials(k), polynomials(pol), ...
                                           points(:, pt_order(point_pt)));
        end
         % Advance pointer to first unused point yet
         point_pt = point_pt + 1;
         if point_pt > n_points
             break
         end
    else
        % No known point satisfies tolerances
        % No op
        1;
    end

end
pivot_polynomials = polynomials;
% Remove unused points
pt_order = pt_order(1:point_pt-1);

end
