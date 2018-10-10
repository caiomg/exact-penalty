function [model, model_changed] = rebuild_model(model, options)
    
    radius_factor = options.radius_factor;
    pivot_threshold = options.pivot_threshold;

    points_abs = [model.points_abs, model.cached_points];
    fvalues = [model.fvalues, model.cached_fvalues];
    radius = model.radius;
    [dim, p_ini] = size(points_abs);
    
    % Center will be first
    points_abs(:, [1, model.tr_center]) = points_abs(:, [model.tr_center, 1]);
    fvalues(:, [1, model.tr_center]) = fvalues(:, [model.tr_center, 1]);
    
    % Recalculate distances 
    points_shifted = zeros(dim, p_ini);
    distances = zeros(1, p_ini);
    for n = 2:p_ini
        points_shifted(:, n) = points_abs(:, n) - points_abs(:, 1);
        distances(n) = norm(points_shifted(:, n), inf); % or 2 norm
    end
    % Reorder
    [distances, pt_order] = sort(distances);
    points_shifted = points_shifted(:, pt_order);
    points_abs = points_abs(:, pt_order);
    fvalues = fvalues(:, pt_order);
    
    % Building model
    pivot_polynomials = band_prioritizing_basis(dim);
    %pivot_polynomials = natural_basis(dim);
    polynomials_num = length(pivot_polynomials);
    pivot_absvalues = zeros(1, polynomials_num);
    % Constant term
    last_pt_included = 1;
    pivot_absvalues(1) = 1;
    poly_i = 2;
    for iter = 2:polynomials_num
        pivot_polynomials(poly_i) = ...
            orthogonalize_to_other_polynomials(pivot_polynomials, ...
                                                    poly_i, ...
                                                    points_shifted, ...
                                                    last_pt_included);
        if poly_i <= dim+1
            % Linear block -- we allow more points
            maxlayer = min(2*radius_factor, distances(end)/radius);
            block_beginning = 2;
            block_end = dim+1;
            if iter > dim+1
                % We already tested all linear terms
                % We do not have points to build a FL model
                break
            end
        else
            % Quadratic block -- being more carefull
            maxlayer = min(radius_factor, distances(end)/radius);
            block_beginning = dim+2;
            block_end = polynomials_num;
        end
        all_layers = linspace(1, maxlayer, ceil(maxlayer));
        max_absval = -inf;
        pt_max = 0;
        for layer = all_layers
            for n = last_pt_included+1:p_ini
                if distances(n) > layer*radius
                    break
                end
                absval = abs(evaluate_polynomial(pivot_polynomials(poly_i), ...
                                                 points_shifted(:, n)));
                if max_absval < absval
                    max_absval = absval;
                    pt_max = n;
                end
            end
        end
        if max_absval > pivot_threshold
            % Point accepted
            pt_next = last_pt_included + 1;
            points_shifted(:, [pt_next, pt_max]) = points_shifted(:, ...
                                                              [pt_max, pt_next]);
            points_abs(:, [pt_next, pt_max]) = points_abs(:, [pt_max, pt_next]);
            fvalues(:, [pt_next, pt_max]) = fvalues(:, [pt_max, ...
                                pt_next]);
            distances([pt_next, pt_max]) = distances([pt_max, pt_next]);
            pivot_absvalues(pt_next) = max_absval;
            
            % Normalize polynomial value
            pivot_polynomials(poly_i) = ...
                normalize_polynomial(pivot_polynomials(poly_i), ...
                                     points_shifted(:, pt_next));
            % Re-orthogonalize
            pivot_polynomials(poly_i) = ...
                orthogonalize_to_other_polynomials(pivot_polynomials, ...
                                                        poly_i, ...
                                                        points_shifted, ...
                                                        last_pt_included);
            
            % Orthogonalize polynomials on present block (deffering
            % subsequent ones)
            pivot_polynomials = orthogonalize_block(pivot_polynomials, ...
                                                    points_shifted(:, ...
                                                              poly_i), ...
                                                    poly_i, ...
                                                    block_beginning, poly_i);
            already_included(pt_max) = true;
            last_pt_included = pt_next;
            poly_i = poly_i + 1;
        else
            % Exchange some polynomials and try to advance...
            % Moving this polynomial to the end of the block
            pivot_polynomials(poly_i:block_end) = ...
                pivot_polynomials([poly_i+1:block_end, poly_i]);
            % If we are on the linear block, this means we won't be
            % able to build a FL model
        end
    
    end
    
    model.tr_center = 1;
    model.points_abs = points_abs(:, 1:last_pt_included);
    model.points_shifted = points_shifted(:, 1:last_pt_included);
    model.fvalues = fvalues(:, 1:last_pt_included);
    model.cached_points = points_abs(:, last_pt_included+1:end);
    model.cached_fvalues = fvalues(:, last_pt_included+1:end);
    model.pivot_polynomials = pivot_polynomials;
    model.pivot_absvalues = pivot_absvalues;
    model.modeling_polynomials = {};
    model_changed = last_pt_included < p_ini;
end


