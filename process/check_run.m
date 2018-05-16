tol_c = 1e-3;
tol_g = 1e-3;
tries = length(all_results);
n_problems = 134
succeeded = false(n_problems, tries);
all_results2 = {};
for k = 1:tries
    this_results = all_results{k};
    succeeded(:, k) = check_run_single(this_results);


    for n = 1:n_problems
        problem = all_problems{n};
        terminate_cutest_problem();
        prob = setup_cutest_problem(problem, '../my_problems/');

        % Objective
        f_obj = @(x) get_cutest_objective(x);

        % Constraints
        n_constraints = get_cutest_total_number_of_constraints();
        all_con = cell(n_constraints, 1);
        for m = 1:n_constraints
            gm = @(x) evaluate_my_cutest_constraint(x, m, 1);
            all_con{m} = gm;
        end
        nlcon = @(x) constraints(all_con, {}, x, 1);
        if succeeded(n, k)
            [is_kkt, lgrad] = check_kkt(f_obj, all_con, this_results(n).x, tol_c, tol_g);
        else
            is_kkt = false;
            lgrad = inf;
        end
        this_results(n).kkt = is_kkt;
        this_results(n).lgrad = lgrad;
        terminate_cutest_problem();
    end
    all_results2{end+1} = this_results;
end
b = sum(succeeded, 2);
bi = [-max(0, min(1, b));
      zeros(tries, 1)];
c = ones(tries, 1);
intcon = 1:tries;
A = [-succeeded;
     -eye(tries)];

chosen = intlinprog(c, intcon, A, bi)