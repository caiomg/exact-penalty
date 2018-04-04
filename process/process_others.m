%results_cobyla = evaluate_real_cutest_functions(results_cobyla);

tol_c = 1e-2;
tol_f = 1e-1;

problems_solved_cobyla = zeros(9000, 1);
problems_solved_l1 = zeros(9000, 12);

for k = 1:n_problems
    sol_fval = results_cobyla(k).sol_fval;
    if isempty(sol_fval)
       sol_fval = results_dfo(k).sol_fval; 
    end
    
    if norm(max(0, results_cobyla(k).real_c)) < tol_c
        if isempty(sol_fval) || results_cobyla(k).fval - sol_fval < tol_f
            n_evals_cobyla = results_cobyla(k).f_count;
            problems_solved_cobyla(n_evals_cobyla:end, 1) = 1 + ...
                problems_solved_cobyla(n_evals_cobyla:end, 1);
        end
    end
    
    
    for m = 1:10
        this_results = df_results(m).results;
        try
            if this_results(k).nphi < tol_c
                if isempty(sol_fval) || this_results(k).fval - sol_fval < tol_f
                    n_evals_this = this_results(k).fcount;
                    problems_solved_l1(n_evals_this:end, m) = 1 + ...
                        problems_solved_l1(n_evals_this:end, m);
                end
            end
        catch
           1; 
        end
    end

        for m = 11:12
        this_results = df_results(m).results;
        try
            if this_results(k).nphi < tol_c
                if isempty(sol_fval) || this_results(k).fx - sol_fval < tol_f
                n_evals_this = this_results(k).fcount;
                problems_solved_l1(n_evals_this:end, m) = 1 + ...
                    problems_solved_l1(n_evals_this:end, m);
                end
            end
        catch
           1; 
        end
    end

end

plot(t, problems_solved_l1);
hold on
plot(t, problems_solved_cobyla, 'k', 'LineWidth', 2);
hold off