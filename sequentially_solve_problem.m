function [these_results, iterations] = sequentially_solve_problem(problem_name, ...
                                                     solver_configuration, ...
                                                     all_mu)


    epsilon = solver_configuration.epsilon;
    Lambda = solver_configuration.Lambda;
    l1_options = solver_configuration.l1_options;
    max_fcount = solver_configuration.max_fcount;
    log_dir = solver_configuration.log_dir;
    feasibility_tol = solver_configuration.feasibility_tol;
    
    [old_iface, prob_iface] = setup_cutest_problem(problem_name, '~/code/problems/');
    dim = prob_iface.dim;
    n_constraints = sum(prob_iface.con_lb > -1e19) + sum(prob_iface.con_ub < 1e19);
    cutest_lower_bounds = prob_iface.lb > -1e19;
    cutest_upper_bounds = prob_iface.con_ub < 1e19;
    

    bl = prob_iface.lb;
    bu = prob_iface.ub;

    fixed_variables = (bl == bu);
    if sum(fixed_variables) == 0
      x0 = prob_iface.x0;
      %%
      lower_bounds = bl > -1e19;
      upper_bounds = bu < 1e19;
      old_constraints_approach = false;
      f_obj = @(x) prob_iface.evaluate_objective(x);
      con_lb = prob_iface.con_lb;
      con_ub = prob_iface.con_ub;
      all_con = cell(prob_iface.n_constraints, 1);
      for q = 1:prob_iface.n_constraints
        all_con{q} = @(x) prob_iface.evaluate_constraint(x, q);
      end 
      else
      lower_bounds = bl > -1e19;
      upper_bounds = bu < 1e19;
      old_constraints_approach = false;
      f_obj = @(x) function_fix_fixed(@(y) prob_iface.evaluate_objective(y), x, bl, bu);
      con_lb = prob_iface.con_lb;
      con_ub = prob_iface.con_ub;
      all_con = cell(prob_iface.n_constraints, 1);
      for q = 1:prob_iface.n_constraints
        all_con{q} = @(x) function_fix_fixed(@(y) prob_iface.evaluate_constraint(y, q), x, bl, bu);
      end 
      bl = bl(~fixed_variables);
      bu = bu(~fixed_variables);
      x0 = prob_iface.x0(~fixed_variables);
    end      
    



    counter = evaluation_counter(f_obj);

    nlcon = @(x) constraints(all_con, {}, x, 1);
    
    % Initial point
    prob_iface.x0;


    f = @(x) counter.evaluate(x);
    these_results = {};
    iterations.name = problem_name;
    iterations.fx = [];
    iterations.violation = [];
    iterations.evaluations = [];
    iterations.radius = [];
    iterations.last_run = 0;
    total_evaluations = 0;
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
            nphi = norm(max(0, max(con_lb - nlcon(x), nlcon(x) - con_ub)));
            [kkt, lgrad] = check_kkt(f, all_con, x, con_lb, con_ub, bl, bu, 1e-5, 5e-5);
            the_error = [];
            iterations.fx = [iterations.fx, hs2.fx];
            iterations.violation = [iterations.violation, hs2.l2_violation];
            iterations.evaluations = [iterations.evaluations, total_evaluations + [hs2.evaluations]];
            iterations.radius = [iterations.radius, hs2.radius];
            iterations.last_run = numel(hs2);
            total_evaluations = iterations.evaluations(end);
        else
            x = [];
            hs2 = [];
            fx = [];
            nphi = [];
            error_x = [];
            fcount = nan;
            kkt = false;
            lgrad = nan;
        end

        problem_result = struct('name', problem_name, ...
                                'x', x, ...
                                'fx', fx, ...
                                'history', [], ...
                                'fcount', fcount, ...
                                'nphi', nphi, ...
                                'mu', mu, ...
                                'epsilon', epsilon, ...
                                'Lambda', Lambda, ...
                                'kkt', kkt, ...
                                'lgrad', lgrad, ...
                                'l1_options', l1_options, ...
                                'the_error', the_error);
        if solved && nphi <= feasibility_tol
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

function [f, grad] = function_fix_fixed(fun, x, lb, ub)
  fixed_variables = (lb == ub);
  x_full = zeros(size(lb));
  x_full(fixed_variables) = lb(fixed_variables);
  x_full(~fixed_variables) = x;
  if nargout < 2
      f = fun(x_full);
  else
    [f, grad] = fun(x_full);
    grad = grad(~fixed_variables);
  end
end

  
