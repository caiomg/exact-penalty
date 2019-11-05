function print_my_table(results)

    n_results = numel(results);
    for k = 1:n_results
        if iscell(results)
            this_result = results{k};
        else
            this_result = results(k);
        end
        print_my_table_row(this_result);
    end

end
