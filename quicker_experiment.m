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

logdir = fullfile('..', 'logs', datestr(now, 30));
if ~exist(logdir, 'dir')
   if mkdir(logdir) ~= 1
       logdir = '.';
   end
end
log_filename = fullfile(logdir, sprintf('%s_p1_db.log', datestr(now, 30)));
log_fd = fopen(log_filename, 'w');

l1_options = [];
l1_options.eta_2 = 0.01;
l1_options.pivot_threshold = 0.01;
% l1_options.basis = 'FULL'
% l1_options.criticality_mu = 0.1 % testing !!!!!
% l1_options.criticality_beta = 0.02 % testing
l1_options.max_iter = 15000;

% Parameters
epsilon = 1;
Lambda = 1e-3;

selected_problems = hs_cheti_sampaio();

all_epsilon = [1]
all_lambda = [0.1]

final_filenames = {};
all_results = {};
good_results = {};
all_solved = [];
n_problems = length(selected_problems);
solved_problems = false(n_problems, 1);
clear results;

solver_configuration.epsilon = epsilon;
solver_configuration.Lambda = Lambda;
solver_configuration.l1_options = l1_options;
solver_configuration.max_fcount = 25000;
solver_configuration.log_dir = logdir;
results = [];

all_mu = 10.^(-1:6);

parfor k = 1:length(selected_problems)

    bad_cond_warn = warning('off', 'cmg:ill_conditioned_system');
    neg_mult_warn = warning('off', 'cmg:multipliers_negative');
    high_mult_warn = warning('off', 'cmg:multipliers_high');
    different_objective_warn = warning('off', 'cmg:different_objective');

    problem_result = handle_problem(selected_problems(k), solver_configuration);

    %     print_result(problem_result, log_fd);
   print_result(problem_result);

   all_results{k} = problem_result;
    
    warning(bad_cond_warn);
    warning(neg_mult_warn);
    warning(high_mult_warn);
    warning(different_objective_warn);
end
for k = 1:length(selected_problems)
    if all_results{k}.kkt || ...
            (~isempty(all_results{k}.nphi) && all_results{k}.nphi < 1e-6 ...
             && (-(all_results{k}.error_rel) < 1e-6 ...
                 ||-(all_results{k}.error_obj) < 1e-7))
        solved_problems(k) = true;
        good_results{end+1} = all_results{k};
        all_solved(end+1) = k;
    end
end

warning('on', 'cmg:badly_conditioned_system');
[~, results_order] = sort(all_solved);
good_results_ordered = {good_results{results_order}};
    filename = fullfile(logdir, sprintf('%s_p1_db', datestr(now, 30)));
    save(filename, 'all_results');
     
print_my_table(good_results_ordered);
problems_not_solved = selected_problems(~solved_problems);
fprintf(1, 'Not solved:\n');
fprintf(1, '  %s ', problems_not_solved.name);
fprintf(1, '\n');