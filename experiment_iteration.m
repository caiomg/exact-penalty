function result = experiment_iteration(problem_name, mu, solution)
    p_seed = rng('default');
    try
        [x, fval, fcount, problem] =  l1_solve_scaled_by_name(problem_name, mu);
        funcs = problem_functions(problem);
        f = funcs{1};
        all_con = {funcs{2:end}};
        nlcon = @(x) constraints(all_con, {}, x, 1);
        con_lb = problem.con_lb;
        con_ub = problem.con_ub;
        lb = problem.lb;
        ub = problem.ub;
        exception_found = [];
        solved = true;
    catch thiserror
        exception_found = thiserror;
        solved = false;
    end
    rng(p_seed);

    if solved
        fx = fval;
        nphi = norm(max(0,max(con_lb - nlcon(x), nlcon(x) - con_ub)));
        error_obj = solution - fx;
        error_rel = error_obj/solution;
        [kkt, lgrad] = check_kkt(f, all_con, x, con_lb, con_ub, lb, ub, 1e-5, 5e-5);
        terminate_cutest_problem(problem.directory);
    else
        x = [];
        hs2 = [];
        fx = [];
        nphi = inf;
        error_obj = inf;
        error_rel = inf;
        error_x = inf;
        fcount = nan;
        kkt = false;
        lgrad = nan;
    end

    result.name = problem_name;
    result.x = x;
    result.fx = fx;
    result.history = []; % hs2
    result.fcount = fcount;
    result.error_obj = error_obj;
    result.error_rel = error_rel;
    result.nphi = nphi;
    result.mu = mu;
%     result.epsilon = epsilon;
%     result.lambda = Lambda;
    result.kkt = kkt;
    result.lgrad = lgrad;
    result.exception_found = exception_found;

end
