function summary = sum_all_results(all_results)

    summary = {};
    n_problems = numel(all_results);
    for k = 1:n_problems
        problem_results = all_results{k};
        total_fcount = 0;
        minimum_violation = inf;
        less_violated = 0;
        for m = 1:numel(problem_results)
            if isfinite(problem_results{m}.fcount)
                total_fcount = total_fcount + problem_results{m}.fcount;
                if problem_results{m}.nphi < minimum_violation
                    minimum_violation = problem_results{m}.nphi;
                    less_violated = m;
                end
            end
        end
        if problem_results{m}.good
            best_result = m;
        else
            best_result = less_violated;
        end
        summary{end+1, 1} = problem_results{best_result};
        summary{end}.best_fcount = summary{end}.fcount;
        summary{end}.fcount = total_fcount;
    end

end