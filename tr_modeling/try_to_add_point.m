function [model, exitflag] = try_to_add_point(model, new_point, new_fvalues, funcs, bl, bu, options)

    point_added = false;
    point_replaced = false;
    
    if ~is_complete(model)
        % Add this point
        relative_pivot_threshold = options.add_threshold;
        [model, point_added] = add_point(model, new_point, new_fvalues, relative_pivot_threshold);
    else
        % Improve model by other point
        [model, point_replaced] = choose_and_replace_point(model, funcs, bl, bu, options);
    end
    if ~point_added && ~point_replaced
        [model, exitflag] = ensure_improvement(model, funcs, bl, bu, options);
    else
        exitflag = 1;
        if isempty(model.modeling_polynomials)
            basis = band_prioritizing_basis(size(model.points_shifted, 1));
            model.modeling_polynomials = ...
                recompute_polynomial_models(model.points_shifted, model.fvalues, ...
                                            basis);
            [~, wid] = lastwarn();
            if strcmp(wid, 'cmg:badly_conditioned_system')
                check = true;
            else
                check = false;
            end
            [dim, points_num] = size(model.points_abs);
            if points_num > dim + 1
               mp1 = shift_polynomial(model.modeling_polynomials{1}, -model.points_abs(:, 1));
               p2 = compute_quadratic_mn_polynomials(model.points_abs, model.tr_center, model.fvalues);
               mp2 = shift_polynomial(p2{1}, -model.points_abs(:, model.tr_center));
               error1 = 0;
               error2 = 0;
               for m = 1:points_num
                   error1 = error1 + abs(evaluate_polynomial(mp1, model.points_abs(:, m)) - model.fvalues(1, m));
                   error2 = error2 + abs(evaluate_polynomial(mp2, model.points_abs(:, m)) - model.fvalues(1, m));
               end
               if check
                   1;
               end
            end
        end
    end
    % Check if model changed

    
end