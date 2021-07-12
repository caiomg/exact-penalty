function print_problems_table(selected_problems, chetis, sampaios)

    for k = 1:numel(selected_problems)
        name = selected_problems(k).name;
        [dim, v_bounds, n_constraints, bound_feasible, nl_feasible] = ...
                                                get_problem_data(name);
        if bound_feasible
            is_feasible_bound = 'yes';
        else
            is_feasible_bound = 'no';
        end
        if nl_feasible
            is_nl_feasible = 'yes';
        else
            is_nl_feasible = 'no';
        end
        in_cheti = strincell(chetis, name);
        in_sampaio = strincell(sampaios, name);
        if in_cheti && in_sampaio
            group = 'A,B';
        elseif in_cheti
            group = 'A';
        elseif in_sampaio
            group = 'B';
        else
            error('cmg:runtime', 'Not in groups');
        end
        
        fprintf(1, ...
            '% 8s  & % 3d  & % 3d  & % 3d  & % 3s  & % 3s  & %s\\\\\n', ...
            name, dim, v_bounds, n_constraints, is_feasible_bound, ...
            is_nl_feasible, group);
            
        
    end


end