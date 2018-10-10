function compute_polynomial_models(points_abs, center_i, fvalues, basis)
    
    
    [dim, points_num] = size(points_abs);
    functions_num = size(fvalues, 1);
    
    linear_terms = dim+1;
    
    if points_num < linear_terms
        % Compute linear model
        
    else
        % Compute quadratic model
        
    end
end
