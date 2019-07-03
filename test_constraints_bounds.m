for k = 1:n_problems
    iter = iter + 1;
%     if solved_problems(k)
%         continue
%     end
    problem_name = selected_problems(k).name;
    prob = setup_cutest_problem(problem_name, '../my_problems/');
    nonlin_min = min(prob.cu - prob.cl);
    bounds_min = min(prob.bu - prob.bl);
    if nonlin_min < 1e10 && bounds_min < 1e10
        problem_name
        nonlin_min
    end
    terminate_cutest_problem();
end