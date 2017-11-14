function all_cutest_constraints = setup_cutest_constraints()
    global problem_data_cutest
    n_constraints = problem_data_cutest.m;
    n_variables = problem_data_cutest.n;
    
    all_cutest_constraints = [];
    cases = 0;
    for k = 1:n_constraints
        if problem_data_cutest.cu(k) < 1e19
            cases = cases + 1;
            all_cutest_constraints(cases).fun = @(x) get_cutest_constraint(x, k);
            all_cutest_constraints(cases).bound = problem_data_cutest.cu(k);
            all_cutest_constraints(cases).sign = 1;
        end
        if problem_data_cutest.cl(k) > -1e19
            cases = cases + 1;
            all_cutest_constraints(cases).fun = @(x) get_cutest_constraint(x, k);
            all_cutest_constraints(cases).bound = problem_data_cutest.cl(k);
            all_cutest_constraints(cases).sign = -1;
        end
    end

    H = zeros(n_variables);
    for k = 1:n_variables
        if problem_data_cutest.bl(k) > -1e19
            cases = cases + 1;
            gk = zeros(n_variables, 1);
            gk(k) = 1;
            this_con = @(x) quadratic(H, gk, 0, x);
            all_cutest_constraints(cases).fun = this_con;
            all_cutest_constraints(cases).bound = problem_data_cutest.bl(k);
            all_cutest_constraints(cases).sign = -1;
        end
        if problem_data_cutest.bu(k) < 1e19
            cases = cases + 1;
            gk = zeros(n_variables, 1);
            gk(k) = 1;
            this_con = @(x) quadratic(H, gk, 0, x);
            all_cutest_constraints(cases).fun = this_con;
            all_cutest_constraints(cases).bound = problem_data_cutest.bu(k);
            all_cutest_constraints(cases).sign = 1;
        end
    end
end