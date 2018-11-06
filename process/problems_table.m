
for k = 1:length(selected_problems)
    problem_name = selected_problems(k).name;
    terminate_cutest_problem(fullfile('../my_problems/', problem_name));
    prob = setup_cutest_problem(problem_name, '../my_problems/');
    dim = prob.n;
    n_constraints = sum(prob.cl > -1e19) + sum(prob.cu < 1e19);
    linear_constraints = sum(prob.linear & prob.cl > -1e19) + sum(prob.linear & prob.cu < 1e19);
    nonlinear_constraints = n_constraints - linear_constraints;
    cutest_lower_bounds = prob.bl > -1e19;
    cutest_upper_bounds = prob.bu < 1e19;
    variable_bounds = sum(cutest_lower_bounds + cutest_upper_bounds);
    terminate_cutest_problem(fullfile('../my_problems/', problem_name));
    fprintf(1, '   % 5s & % 3d   & % 3d   & % 3d   & % 3d \\\\\n', problem_name, dim, nonlinear_constraints, linear_constraints, variable_bounds);
end