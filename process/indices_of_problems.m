function ind = indices_of_problems(results, problem_names)

    n_names = numel(problem_names);
    n_results = numel(results);
    result_names = {results.name};
    for n = 1:n_names
        current_name = problem_names{n};
        name_found = false;
        for r = 1:n_results
            if strcmp(result_names{r}, current_name)
                name_found = true;
                break
            end
        end
        if name_found
            ind(n) = r;
        else
            % either raise error or an invalid index
            %error('cmg:notfound', 'Name not found');
            ind(n) = nan;
        end
    end

end