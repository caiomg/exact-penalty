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
if exist('l1_solve_scaled_by_name', 'file') ~= 2
   addpath('../scaling_problems'); 
end

selected_problems = generate_selected_problems();

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
l1_options.basis = 'FULL'
%l1_options.crit_mu = 0.1 % testing !!!!!

% Parameters
epsilon = 1;
delta = 1e-6;
Lambda = 0.1;


all_epsilon = [1]
all_lambda = [0.1]


final_filenames = {};
all_results = {};
good_results = [];
all_solved = [];
n_problems = length(selected_problems);
problem_indices = 1:n_problems;
solved_problems = false(n_problems, 1);
mu_used = inf(n_problems, 1);
result_pt = 0;
clear results;
results = struct('name', {}, ...
                 'x', {}, ...
                 'fx', {}, ...
                 'history', {}, ...
                 'fcount', {}, ...
                 'error_obj', {}, ...
                 'error_rel', {}, ...
                 'nphi', {}, ...
                 'mu', {}, ...
                 'kkt', {}, ...
                 'lgrad', {}, ...
                 'exception_found', {});


all_mu = [1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6];

tic
for mu = all_mu
    initially_not_solved = ~solved_problems;
    problems_to_solve = problem_indices(initially_not_solved);
    n_not_solved = numel(problems_to_solve);
    last_index = numel(results);
    pass_status = false(n_not_solved, 1);
    parfor k = 1:n_not_solved
        bad_cond_warn = warning('off', 'cmg:ill_conditioned_system');
        fumou_warn =  warning('off', 'cmg:fumou');
        ascent_warn = warning('off', 'cmg:not_descent');
        prob_index = problems_to_solve(k);
        problem_name = selected_problems(prob_index).name;
        result = experiment_iteration(problem_name, mu, selected_problems(prob_index).solution);

        if result.kkt || (result.nphi < 1e-6 && -result.error_rel < 1e-6)
            pass_status(k) = true;
        end
        results(last_index + k) = result;

%         print_results(result, log_fd);
        print_results(result, 1);
    end
    for k = 1:n_not_solved
        print_results(results(last_index + k), log_fd);
    end
    problems_solved_in_pass = problems_to_solve(pass_status);
    solved_in_pass = false(n_problems, 1);
    solved_in_pass(problems_solved_in_pass) = true(1, sum(pass_status));
    all_solved = [all_solved, problems_solved_in_pass];
    successful_indices = 1:n_not_solved;
    successful_indices = last_index + successful_indices(pass_status);
    good_results = [good_results, results(successful_indices)]; 
    solved_problems = solved_problems | solved_in_pass;
    fprintf(1, ['----------------------------------------', ...
                '----------------------------------------', '\n']);
end
toc
[~, results_order] = sort(all_solved);
good_results_ordered = good_results(results_order);
    filename = fullfile(logdir, sprintf('%s_p1_db', datestr(now, 30)));
    save(filename, 'results', 'mu', 'epsilon', ...
         'delta', 'Lambda', 'final_filenames', 'good_results', ...
         'all_solved', 'good_results_ordered');
     
print_my_table(good_results_ordered);
