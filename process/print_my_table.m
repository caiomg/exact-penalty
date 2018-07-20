function print_my_table(results)

    n_results = length(results);
    for k = 1:n_results
        this_result = results{k};
        print_my_table_row(this_result);
    end
    