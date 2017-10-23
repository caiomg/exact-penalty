function m = get_cutest_total_number_of_constraints()
    global problem_data_cutest
    n_constraints = problem_data_cutest.m;
    n_variables = problem_data_cutest.n;

    % In the presence of upper bound constraints
    if min(problem_data_cutest.bu) < 1e15
        m = n_constraints + n_variables + n_variables;
    elseif max(problem_data_cutest.bl) > -1e15
        m = n_constraints + n_variables;
    else
        m = n_constraints;
    end
end