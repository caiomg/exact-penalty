function [new_point, new_pivot_value, point_found] = ...
        point_new(polynomial, x_tr_center, radius, bl, bu, pivot_threshold)
    
    
    [new_point_min, pivot_min, exitflag_min] = minimize_tr(polynomial, ...
                                                      x_tr_center, ...
                                                      radius, bl, bu);
    
    polynomial_max = multiply_p(polynomial, -1);
    [new_point_max, pivot_max, exitflag_max] = minimize_tr(polynomial_max, ...
                                                      x_tr_center, ...
                                                      radius, bl, ...
                                                      bu);
    if (exitflag_min == 0 || exitflag_min == 1) && (abs(pivot_min) ...
                                                    > pivot_threshold)
        if (exitflag_max == 0 || exitflag_max == 1) && (abs(pivot_max) ...
                                                        >= abs(pivot_min))
            new_point = new_point_max;
            new_pivot_value = pivot_max;
            point_found = true;
        else
            new_point = new_point_min;
            new_pivot_value = pivot_min;
            point_found = true;
        end
    elseif (exitflag_max == 0 || exitflag_max == 1) && (abs(pivot_max) ...
                                                        >= pivot_threshold)
        new_point = new_point_max;
        new_pivot_value = pivot_max;
        point_found = true;
    else
        point_found = false;
        new_point = [];
        new_pivot_value = 0;
    end
    
        
end
