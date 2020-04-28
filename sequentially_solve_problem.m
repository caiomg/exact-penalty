function these_results = sequentially_solve_problem(problem_data, ...
                                                     solver_configuration, ...
                                                     all_mu)
% HANDLE_PROBLEM - 
%   


    problem_name = problem_data.name;
    epsilon = solver_configuration.epsilon;
    Lambda = solver_configuration.Lambda;
    l1_options = solver_configuration.l1_options;
    max_fcount = solver_configuration.max_fcount;
    log_dir = solver_configuration.log_dir;
    
    [old_iface, prob_iface] = setup_cutest_problem(problem_name, '../my_problems');
    dim = prob_iface.dim;
    n_constraints = sum(prob_iface.con_lb > -1e19) + sum(prob_iface.con_ub < 1e19);
    cutest_lower_bounds = prob_iface.lb > -1e19;
    cutest_upper_bounds = prob_iface.con_ub < 1e19;
    

    bl = prob_iface.lb;
    bu = prob_iface.ub;
%%
    lower_bounds = bl > -1e19;
    upper_bounds = bu < 1e19;
    old_constraints_approach = false;
    if old_constraints_approach
        % Objective
        f_obj = @(x) get_cutest_objective(x);
        % 'Nonlinear' constraints
        all_con = cell(n_constraints, 1);
        for q = 1:n_constraints
            gk = @(x) evaluate_my_cutest_constraint(x, q, 1);
            all_con{q} = gk;
        end
        len_con = length(all_con);
        con_lb = -inf(len_con, 1);
        con_ub = zeros(len_con, 1);
    else
        f_obj = @(x) prob_iface.evaluate_objective(x);
        con_lb = prob_iface.con_lb;
        con_ub = prob_iface.con_ub;
        con_lb(con_lb < -1e19) = -inf(size(con_lb(con_lb < -1e19)));
        con_ub(con_ub > 1e19) = inf(size(con_ub(con_ub > 1e19)));
        all_con = cell(prob_iface.n_constraints, 1);
        for q = 1:prob_iface.n_constraints
            all_con{q} = @(x) prob_iface.evaluate_constraint(x, q);
        end 
    end
    counter = evaluation_counter(f_obj);

    nlcon = @(x) constraints(all_con, {}, x, 1);
    
    % Initial point
    prob_iface.x0;


    f = @(x) counter.evaluate(x);
    x0 = prob_iface.x0;
    these_results = {};
    for mu = all_mu
        counter.reset_count();
        counter.set_max_count(max_fcount);

        try
            p_seed = rng('default');

            [x, hs2] = l1_penalty_solve(f, all_con, con_lb, con_ub, x0, ...
                                        mu, epsilon, [], Lambda, bl, ...
                                        bu, l1_options);
            solved = true;
        catch thiserror
            the_error = thiserror;
            solved = false;
        end
        rng(p_seed);
        if solved
            fcount = counter.get_count();
            fx = f(x);
            nphi = norm(max(0,max(con_lb - nlcon(x), nlcon(x) - con_ub)));
            error_obj = problem_data.solution - fx;
            error_rel = -error_obj/(hs2(1).fx - problem_data.solution);
            [kkt, lgrad] = check_kkt(f, all_con, x, con_lb, con_ub, bl, bu, 1e-5, 5e-5);
            the_error = [];
        else
            x = [];
            hs2 = [];
            fx = [];
            nphi = [];
            error_obj = [];
            error_x = [];
            fcount = nan;
            kkt = false;
            error_rel = inf;
            lgrad = nan;
        end

        problem_result = struct('name', problem_name, ...
                                'x', x, ...
                                'fx', fx, ...
                                'history', [], ...
                                'fcount', fcount, ...
                                'error_obj', error_obj, ...
                                'nphi', nphi, ...
                                'mu', mu, ...
                                'epsilon', epsilon, ...
                                'Lambda', Lambda, ...
                                'kkt', kkt, ...
                                'error_rel', error_rel, ...
                                'lgrad', lgrad, ...
                                'l1_options', l1_options, ...
                                'the_error', the_error);
        if problem_result.kkt ...
                || (~isempty(problem_result.nphi) && (problem_result.nphi < 1e-5) ...
                    && (-problem_result.error_obj <  1e-6 ...
                        || problem_result.error_rel < 1e-4))
            problem_result.good =  true;
        else
            problem_result.good = false;
        end

        filename = fullfile(log_dir, [problem_name, '.result']);
        fid = fopen(filename, 'a');
        fprintf(fid, '\nmu:%d     kkt: %d', mu, kkt); 
        fprintf(fid, '\n%s\n', jsonencode(problem_result));
        fclose(fid);
        these_results{1,end+1} = problem_result;
        if problem_result.good
            break
        end

    end
        terminate_cutest_problem();
    
end
