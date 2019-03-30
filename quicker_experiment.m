if exist('terminate_cutest_problem', 'file') ~= 2
  addpath('my_cutest_functions');
end
if exist('check_kkt', 'file') ~= 2
  addpath('process');
end
if exist('evaluate_polynomial', 'file') ~= 2
  addpath('tr_modeling');
  addpath('tr_modeling/polynomials');
end



terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest

logdir = fullfile('..', 'logs', datestr(now, 30));
if ~exist(logdir, 'dir')
   if mkdir(logdir) ~= 1
       logdir = '.';
   end
end
log_filename = fullfile(logdir, sprintf('%s_p1_db.log', datestr(now, 30)));
log_fd = fopen(log_filename, 'w');


l1_options = [];
l1_options.eta_2 = 0.05;
l1_options.pivot_threshold = 0.001;
l1_options.basis = 'dummy'

% Parameters
epsilon = 1;
delta = 1e-6;
Lambda = 0.1;

list_of_problems


all_epsilon = [1]
all_lambda = [0.1]


final_filenames = {};
all_results = {};
good_results = {};
all_solved = [];
iter = 0;
n_problems = length(selected_problems);
solved_problems = false(n_problems, 1);
clear results;

bad_cond_warn = warning('off', 'cmg:badly_conditioned_system');

reordered_problems = [1:14, 16:30, 15];

for k = reordered_problems
    iter = iter + 1;

    problem_name = selected_problems(k).name;
    mu = selected_problems(k).mu;
    [prob, prob_iface] = setup_cutest_problem(problem_name, '~/docd/exchange/my_problems');
    dim = prob.n;
    n_constraints = sum(prob.cl > -1e19) + sum(prob.cu < 1e19);
    cutest_lower_bounds = prob.bl > -1e19;
    cutest_upper_bounds = prob.bu < 1e19;
    

    bl = prob.bl;
    bu = prob.bu;
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
        con_lb = prob.cl;
        con_ub = prob.cu;
        con_lb(con_lb < -1e19) = -inf(size(con_lb(con_lb < -1e19)));
        con_ub(con_ub > 1e19) = inf(size(con_ub(con_ub > 1e19)));
        all_con = cell(prob.m, 1);
        for q = 1:prob.m
            all_con{q} = @(x) prob_iface.evaluate_constraint(x, q);
        end 
    end
    counter = evaluation_counter(f_obj);

    nlcon = @(x) constraints(all_con, {}, x, 1);
    
    % Initial point
    x0 = prob.x;


    fixed_scale = (prob.bu - prob.bl)/2;
    no_scale = isinf(fixed_scale);

    no_scale = true(size(x0)); % Removing scale
    fixed_scale(no_scale) = no_scale(no_scale);
    O = (prob.bu + prob.bl)/2;
    O(no_scale) = x0(no_scale);
    O = zeros(size(x0)); % Removing scale
    s0 = (x0 - O)./fixed_scale;
    Sc = diag(fixed_scale);
    
    for q = 1:length(all_con)
       all_con{q} = @(w) scale_function(@(x) all_con{q}(x), O, Sc, w);
    end
%     f = @(x) counter.evaluate(x);
    f = @(w) scale_function(@(x) counter.evaluate(x), O, Sc, w);

    


    counter.reset_count();
    counter.set_max_count(50000);

    try
        p_seed = rng('default');

        [sf, hs2] = l1_penalty_solve(f, all_con, con_lb, con_ub, s0, ...
                                     mu, epsilon, delta, Lambda, bl, ...
                                     bu, l1_options);
        x = O + Sc*sf;
        solved = true;
    catch thiserror
        results(k, 1).except = thiserror;
        solved = false;
    end
    rng(p_seed);
    if solved
        fcount = counter.get_count();
        fx = f(x);
        nphi = norm(max(0,max(con_lb - nlcon(x), nlcon(x) - con_ub)));
        error_obj = selected_problems(k).solution - fx;
        error_rel = error_obj/abs(selected_problems(k).solution);
        [kkt, lgrad] = check_kkt(f, all_con, x, con_lb, con_ub, bl, bu, 1e-5, 5e-5);
    else
        x = [];
        hs2 = [];
        fx = [];
        nphi = [];
        error_obj = [];
        error_x = [];
        fcount = nan;
        kkt = false;
        lgrad = nan;
    end

    results(k, 1).name = problem_name;
    results(k, 1).x = x;
    results(k, 1).fx = fx;
    results(k, 1).history = []; % hs2
    results(k, 1).fcount = fcount;
    results(k, 1).error_obj = error_obj;
    results(k, 1).nphi = nphi;
    results(k, 1).mu = mu;
    results(k, 1).epsilon = epsilon;
    results(k, 1).lambda = Lambda;
    results(k, 1).kkt = kkt;
    results(k, 1).lgrad = lgrad;

    print_results(results(k, 1), log_fd);
    print_results(results(k, 1));

    all_results{iter} = results;
    if kkt || ...
            (~isempty(nphi) && nphi < 1e-6 && -error_rel < 1e-6)
        solved_problems(k) = true;
        good_results{sum(solved_problems)} = results(k, 1);
        all_solved(end+1) = k;
    end
    terminate_cutest_problem();
end
warning('on', 'cmg:badly_conditioned_system');
[~, results_order] = sort(all_solved);
good_results_ordered = {good_results{results_order}};
    filename = fullfile(logdir, sprintf('%s_p1_db', datestr(now, 30)));
    save(filename, 'all_results', 'mu', 'epsilon', ...
         'delta', 'Lambda', 'final_filenames', 'good_results', ...
         'all_solved', 'good_results_ordered');
     
print_my_table(good_results_ordered);
