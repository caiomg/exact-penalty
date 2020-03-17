function computing_times = harder_problems(params, indices)
% HARDER_PROBLEMS - 
%   


    l1_options = [];
    l1_options.basis = 'FULL';
    l1_options.crit_mu = 0.1; % testing !!!!!

    % Parameters
    epsilon = params(1);
    Lambda = params(2);
    l1_options.eta_2 = params(3);
    % l1_options.pivot_threshold = 0.001;
    l1_options.pivot_threshold = params(4);

    selected_problems = hs_cheti_problems();

    all_results = {};

    solver_configuration.epsilon = epsilon;
    solver_configuration.Lambda = Lambda;
    solver_configuration.l1_options = l1_options;
    solver_configuration.max_fcount = 10000;
    results = [];
    selected_problems = selected_problems(indices);
    n_problems = length(selected_problems);
    parfor k = 1:n_problems

        bad_cond_warn = warning('off', 'cmg:ill_conditioned_system');
        neg_mult_warn = warning('off', 'cmg:multipliers_negative');
        high_mult_warn = warning('off', 'cmg:multipliers_high');

        problem_result = handle_problem(selected_problems(k), solver_configuration);

        %     print_result(problem_result, log_fd);
%         print_result(problem_result);

        all_results{k} = problem_result;
        
        warning(bad_cond_warn);
        warning(neg_mult_warn);
        warning(high_mult_warn);
    end

    computing_times = nan(n_problems, 1);
    for k = 1:n_problems
        if all_results{k}.kkt ...
                || (~isempty(all_results{k}.error_obj) ...
                    && all_results{k}.error_obj > -1e-5 ...
                    && all_results{k}.nphi < 1e-5)
            computing_times(k) = all_results{k}.fcount;
        else
            computing_times(k) = nan;
        end
    end
end
