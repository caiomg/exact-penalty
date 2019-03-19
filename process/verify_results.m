function summary = verify_results(summary)

    directory = '~/docd/exchange/my_problems';
    bigM = 1e19;
    for k = 1:length(summary)
        if isempty(summary{k}.x)
            summary{k}.fx = [];
            summary{k}.convals = [];
            summary{k}.lgrad = [];
        else
            prob_name = summary{k}.name;
            prob_path = fullfile(directory, prob_name);
            terminate_cutest_problem(prob_path)

            prob = setup_cutest_problem(prob_name, directory);

            f = @(x) get_cutest_objective(x);
            summary{k}.fx = f(summary{k}.x);

            n_constraints = sum(prob.cl > -bigM) + sum(prob.cu < bigM);
            l_bounds = (prob.bl > -bigM);
            u_bounds = (prob.bu < bigM);
            cvals = zeros(n_constraints, 1);

            all_con = cell(prob.m, 1);
            for q = 1:n_constraints
                all_con{q} = @(x) evaluate_my_cutest_constraint(x, q, 1);
                cvals(q) = all_con{q}(summary{k}.x);
            end
            lbound_vals = prob.bl(l_bounds) - summary{k}.x(l_bounds);
            ubound_vals = summary{k}.x(u_bounds) - prob.bu(u_bounds);

            summary{k}.convals = [cvals; lbound_vals; ubound_vals];

            [kkt_ok, lgrad] = check_kkt(f, all_con, summary{k}.x, prob.bl, prob.bu, 1e-5, 5e-5);
            summary{k}.lgrad = lgrad;

            terminate_cutest_problem(prob_path)
        end
    end



end