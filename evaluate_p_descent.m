function descent = evaluate_p_descent(old_values, new_values, con_lb, ...
                                      con_ub, mu)
    
    f_change = old_values(1) - new_values(1);
    
    old_values = old_values(2:end);
    new_values = new_values(2:end);
    
    old_cvalues = [con_lb(isfinite(con_lb)) - old_values(isfinite(con_lb));
                   old_values(isfinite(con_ub)) - con_ub(isfinite(con_ub))];
                   
    new_cvalues = [con_lb(isfinite(con_lb)) - new_values(isfinite(con_lb));
                   new_values(isfinite(con_ub)) - con_ub(isfinite(con_ub))];
    
    c_change = max(0, old_cvalues) - max(0, new_cvalues);
    
    descent = f_change + mu*sum(sort(c_change));
                   
end