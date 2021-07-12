
unfeasible = struct('name', {}, 'max_violation', {}, 'total_violation', {});

for n = 1:numel(selected_problems)
    problem = selected_problems(n);
    [~, p] = setup_cutest_problem(problem.name, '~/code/my_problems');
    x = project_to_bounds(p.x0, p.lb, p.ub);
    max_violation = 0;
    total_violation = 0;
    for k = 1:p.n_constraints
        constraint_value = p.evaluate_constraint(x, k);
        violation = max([0, p.con_lb(k) - constraint_value, constraint_value - p.con_ub(k)]);
        max_violation = max(violation, max_violation);
        total_violation = total_violation + violation;
    end
    if max_violation > 0
        this_problem.name = problem.name;
        this_problem.max_violation = max_violation;
        this_problem.total_violation = total_violation;
        unfeasible(end+1) = this_problem;
    end
end