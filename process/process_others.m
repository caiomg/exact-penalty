%results_cobyla = evaluate_real_cutest_functions(results_cobyla);

tol_c = 1e-4;
tol_f = 1e-1;

tmax = 9000;
problems_solved_cobyla = zeros(tmax, 1);
problems_solved_l1 = zeros(tmax, 12);

cobyla_solved = false(n_problems, 1);
l1_solved = false(n_problems, 12);
cobyla_better = zeros(n_problems, 12);
cobyla_evals = tmax + zeros(n_problems, 1);
l1_evals = tmax + zeros(n_problems, 12);

for k = 1:n_problems
    sol_fval = [];
    %     if isempty(sol_fval)
%        sol_fval = results_dfo(k).sol_fval; 
%     end
    
    if norm(max(0, results_cobyla(k).real_c)) < tol_c
        n_evals_cobyla = results_cobyla(k).f_count;
        problems_solved_cobyla(n_evals_cobyla:end, 1) = 1 + ...
            problems_solved_cobyla(n_evals_cobyla:end, 1);
        sol_fval = results_cobyla(k).sol_fval;
        cobyla_solved(k) = true;
        cobyla_evals(k) = results_cobyla(k).f_count;
    end
    
    
    for m = 1:12
        this_results = all_results2{m};
        if this_results(k).nphi < tol_c
            n_evals_this = this_results(k).fcount;
            problems_solved_l1(n_evals_this:end, m) = 1 + ...
                problems_solved_l1(n_evals_this:end, m);
            l1_solved(k, m) = true;
            l1_evals(k, m) = this_results(k).fcount;
        end
        if cobyla_solved(k)
            if l1_solved(k, m)
                if results_cobyla(k).sol_fval < this_results(k).fx - tol_f
                    cobyla_better(k, m) = 1;
                elseif this_results(k).fx < results_cobyla(k).sol_fval - tol_f
                    cobyla_better(k, m) = -1;
                end
            else
                cobyla_better(k, m) = 1;
            end
        elseif l1_solved(k, m)
            cobyla_better(k, m) = -1;
        end
            
    end
    
    
% 
%     for m = 11:12
%         this_results = df_results(m).results;
%         try
%             if this_results(k).nphi < tol_c
%                 if isempty(sol_fval) || this_results(k).fx - sol_fval < tol_f
%                 n_evals_this = this_results(k).fcount;
%                 problems_solved_l1(n_evals_this:end, m) = 1 + ...
%                     problems_solved_l1(n_evals_this:end, m);
%                 end
%             end
%         catch
%            1; 
%         end
%     end

end

t = 1:tmax;
plot(t, problems_solved_l1(:, [8, 9, 11, 12]));
hold on
plot(t, problems_solved_cobyla, 'k', 'LineWidth', 2);
hold off

any_solved = sum([cobyla_solved, l1_solved(:, :)], 2) > 0;
[rho, tau] = dm_performance_profile([cobyla_evals, l1_evals], any_solved)