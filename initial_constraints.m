
for k = 1:numel(selected_problems)
    problem_name = selected_problems(k).name;
    [~, prob_iface] = setup_cutest_problem(problem_name, '../my_problems');    
    x0 = prob_iface.x0;
    con_lb = prob_iface.con_lb;
    con_ub = prob_iface.con_ub;
    con_lb(con_lb < -1e19) = -inf(size(con_lb(con_lb < -1e19)));
    con_ub(con_ub > 1e19) = inf(size(con_ub(con_ub > 1e19)));
    n_constraints = prob_iface.n_constraints;
    con_values = nan(2*n_constraints, 1);
    for q = 1:n_constraints
        con_val = prob_iface.evaluate_constraint(x0, q);
        con_values(q) = abs(con_lb(q) - con_val);
        con_values(2*q) = abs(con_ub(q) - con_val);
    end
    positives = con_values > 0;
    if ~isempty(find(positives, 1))
        min_positives(k) = min(con_values(positives));
    else
        min_positives(k) = 0;
    end
    considered = isfinite(con_values);
    mean_cons(k) = mean(con_values(considered));
    median_cons(k) = median(con_values(considered));
    fprintf(1, '| % 7s |  % 10g |  % 10g |  % 10g |\n', ...
            problem_name, min_positives(k), mean_cons(k), median_cons(k));
end


    