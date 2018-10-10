function [points_perm, pivot_absvalues, quadratic_pivot_polynomials, ...
          polynomials_permutation] = check_geometry_quadratic_mn(points, ...
                                                      basis_q, tol_pivot)
    

    [dim, points_num] = size(points);
    % Linear part
    % Each iteration checks a point for the given function
    pivot_absvalues = zeros(1, length(basis_q));
    points_perm = 1:points_num;
    poly_perm = 1:length(basis_q);
    poly_i = 1;
    point_i = 1;
    for k = 1:length(basis_q)
        pol_current = basis_q(poly_perm(poly_i));
        
        % Find good pivot
        max_abs_pt = point_i;
        max_abs = -1;
        for p = point_i:points_num
            point_current = points(:, points_perm(p));
            abs_value = abs(evaluate_polynomial(pol_current, ...
                                                point_current));
            if abs_value > max_abs
                max_abs = abs_value;
                max_abs_pt = p;
            end
        end
        if max_abs > tol_pivot
            % Exchange points
            % Maybe I should only exchange when the current point
            % is not that close...
            points_perm([point_i, max_abs_pt]) = ...
                points_perm([max_abs_pt, point_i]);
            point_current = points(:, points_perm(point_i));
            % Eliminate
            for m = poly_i+1:length(basis_q)
                basis_q(poly_perm(m)) = zero_at_point(basis_q(poly_perm(m)), ...
                                           pol_current, point_current);
            end
            pivot_absvalues(point_i) = max_abs;
            % Advance indices...
            point_i = point_i + 1;
            poly_i = poly_i + 1;
        else
            % Maybe just break...
            poly_perm(poly_i:end) = [poly_perm(poly_i+1:end), poly_perm(poly_i)];
        end
    end
    quadratic_pivot_polynomials = basis_q;
    polynomials_permutation = poly_perm;