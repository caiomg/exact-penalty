function print_my_table(all_results)

    n_results = length(all_results);
    for k = 1:n_results
        this_result = all_results{k}(end);
        print_my_table_row(this_result);
    end
    