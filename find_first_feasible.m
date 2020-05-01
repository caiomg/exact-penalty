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
l1_options.pivot_threshold = 0.04
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
all_feasibility_info = {};
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

    feasibility_info = handle_problem_first_feasible(selected_problems(k), solver_configuration);

    %     print_result(problem_result, log_fd);

   all_feasibility_info{k} = feasibility_info;
    
    warning(bad_cond_warn);
    warning(neg_mult_warn);
    warning(high_mult_warn);
    warning(different_objective_warn);
end

