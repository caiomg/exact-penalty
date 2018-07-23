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


l1_options = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_1', 0, 'eta_2', 0.1, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 0.5, 'radius_max', 1e3, ...
                        'criticality_mu', 50, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'diagonal hessian', ...
                        'pivot_threshold', 0.1, 'poised_radius_factor', 2, ...
                        'pivot_imp', 1.1)

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

fprintf(log_fd, ['\nmu = % 6d,    epsilon = % 8g,    delta = % 8g,    ' ...
                 'Lambda = % 8g\n'], mu, epsilon, delta, Lambda);
fprintf(1, ['\nmu = % 6d,    epsilon = % 8g,    delta = % 8g,    ' ...
                 'Lambda = % 8g\n'], mu, epsilon, delta, Lambda);


for k = 15;%1:n_problems
    iter = iter + 1;
%     if solved_problems(k)
%         continue
%     end
    problem_name = selected_problems(k).name;
    mu = selected_problems(k).mu;
    prob = setup_cutest_problem(problem_name, '../my_problems/');
    dim = prob.n;
    n_constraints = sum(prob.cl > -1e19) + sum(prob.cu < 1e19);
    cutest_lower_bounds = prob.bl > -1e19;
    cutest_upper_bounds = prob.bu < 1e19;
    
    % Objective
    f_obj = @(x) get_cutest_objective(x);
    counter = evaluation_counter(f_obj);
    f = @(x) counter.evaluate(x);
    bl = [];
    bu = [];
    % bl = prob.bl;
    % bu = prob.bu;
%%
    lower_bounds = bl > -1e19;
    upper_bounds = bu < 1e19;
    % 'Nonlinear' constraints
    all_con = cell(n_constraints, 1);

    for q = 1:n_constraints
        gk = @(x) evaluate_my_cutest_constraint(x, q, 1);
        all_con{q} = gk;
    end

    for q = 1:dim
        if prob.bl(q) > -1e19
            % Lower bound to consider
            if (isempty(bl) || bl(q) < -1e19)
                % Include constraint as black-box function
                H = zeros(dim);
                g = zeros(dim, 1);
                g(q) = -1;
                c = prob.bl(q);
                gk = @(x) quadratic(H, g, c, x);
                all_con{end+1} = gk;
            else
                % Considered explicitly
                % pass
                1;
            end
        end
        if prob.bu(q) < 1e19
            % Upper bound to consider
            if (isempty(bu) || bu(q) < -1e19)
                % Include constraint as black-box function
                H = zeros(dim);
                g = zeros(dim, 1);
                g(q) = 1;
                c = -prob.bu(q);
                gk = @(x) quadratic(H, g, c, x);
                all_con{end+1} = gk;
            else
                % Considered explicitly
                % pass
                1;
            end
        end
    end

    if length (all_con) ~= n_constraints + sum(cutest_lower_bounds) + ...
            sum(cutest_upper_bounds) - sum(lower_bounds) - sum(upper_bounds)
        error();
    end

    % bl = prob.bl;
    % bu = prob.bu;

%%

    % Initial point
    x0 = prob.x;

    nlcon = @(x) constraints(all_con, {}, x, 1);
%                 fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
%                                                'SpecifyObjectiveGradient', true);
%                 [x_fmincon, fx_fmincon, exitflag, output, mult_fmincon] = fmincon(f, x0,[],[],[],[],bl,bu, nlcon, fmincon_options);
%                 fprintf(1, '| %8s |', problem_name);
%                 fprintf(1, '  % +9.3g |\n', sum(abs([mult_fmincon.lower; mult_fmincon.upper; mult_fmincon.ineqnonlin])));

%                 fcount_fmincon = counter.get_count();
%                 fx_fmincon = f(x_fmincon);
%                 nphi_fmincon = norm(max(0, nlcon(x_fmincon)));
    %%

    counter.reset_count();
    counter.set_max_count(10000);

    try
        p_seed = rng('default');
        [x, hs2] = l1_penalty_solve(f, all_con, x0, mu, epsilon, delta, Lambda, bl, bu, l1_options);
        solved = true;
    catch thiserror
        results(k, 1).except = thiserror;
        solved = false;
    end
    rng(p_seed);
    if solved
        fcount = counter.get_count();
        fx = f(x);
        nphi = norm(max(0, nlcon(x)));
        error_obj = selected_problems(k).solution - fx;
        [kkt, lgrad] = check_kkt(f, all_con, x, bl, bu, 1e-6, 1e-4);
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
            (~isempty(nphi) && nphi < 1e-6 && abs(error_obj) < 5e-6)
        solved_problems(k) = true;
        good_results{sum(solved_problems)} = results(k, 1);
        all_solved(end+1) = k;
    end
    terminate_cutest_problem();
end
[~, results_order] = sort(all_solved);
good_results_ordered = {good_results{results_order}};
    filename = fullfile(logdir, sprintf('%s_p1_db', datestr(now, 30)));
    save(filename, 'all_results', 'mu', 'epsilon', ...
         'delta', 'Lambda', 'final_filenames', 'good_results', ...
         'all_solved', 'good_results_ordered');
     
print_my_table(good_results_ordered);
