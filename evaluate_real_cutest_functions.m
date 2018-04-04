function results = evaluate_real_cutest_functions(results)

n_problems = length(results);
terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest

for n = 1:n_problems

    problem_name = results(n).name;
    prob = setup_cutest_problem(problem_name, fullfile('..' , 'my_problems'));

    % Objective
    f_obj = @(x) get_cutest_objective(x);

    % Constraints
    n_constraints = get_cutest_total_number_of_constraints();
    all_con = cell(n_constraints, 1);
    for k = 1:n_constraints
        gk = @(x) evaluate_my_cutest_constraint(x, k, 1);
        all_con{k} = gk;
    end

    nlcon = @(x) constraints(all_con, {}, x, 1);

    results(n).real_f = f_obj(results(n).x);
    results(n).real_c = nlcon(results(n).x);
    terminate_cutest_problem()
    
end
