
terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest

problem_name = 'HS100';
% problem_name = 'WOMFLET';
% problem_name = 'POLAK1';
% problem_name = 'SNAKE';
% problem_name = 'POLAK3'; % Second derivative of cons. too much important
% % problem_name = 'QC';
% % problem_name = 'CB2';
% % problem_name = 'LOOTSMA';
% % problem_name = 'HS88';
% problem_name = 'HS18';
% problem_name = 'HS19';
% problem_name = 'HS21';
% problem_name = 'HS101';
problem_name = 'HS98';

prob = setup_cutest_problem(problem_name, '../my_problems/');

% Objective
f_obj = @(x) get_cutest_objective(x);
counter = evaluation_counter(f_obj);
f = @(x) counter.evaluate(x);

% Constraints
n_constraints = get_cutest_total_number_of_constraints();

bl = [];
bu = [];
% Bound constraints
lower_bounds = prob.bl > -1e19;
upper_bounds = prob.bu < 1e19;
bl = prob.bl;
bu = prob.bu;
% Remove constraints that actually are bounds
n_constraints = n_constraints - sum(lower_bounds) - sum(upper_bounds);

% Other constraints
all_con = cell(n_constraints, 1);
for k = 1:n_constraints
    gk = @(x) evaluate_my_cutest_constraint(x, k, 1);
    all_con{k} = gk;
end

% Initial point
x0 = prob.x;

% Parameters
mu = 750;

epsilon = 0.85;
delta = 1e-6;
Lambda = 0.075;


nlcon = @(x) constraints(all_con, {}, x, 1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);

x_fmincon = fmincon(f, x0,[],[],[],[], bl, bu, nlcon, fmincon_options);
fx_fmincon = f(x_fmincon);
nlcon_fmincon = max(0, nlcon(x_fmincon));

counter.get_count()
counter.reset_count()
counter.set_max_count(15000);
%%
p_seed = rng('default');
% [x, hs2] = l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda)
[x, hs2] = l1_penalty_solve(f, all_con, x0, mu, epsilon, delta, Lambda, bl, bu, [])
fx = f(x)
nphi = norm(max(0, nlcon(x)))

rng(p_seed);

counter.get_count()
counter.reset_count()


tl1 = @() l1_penalty_solve(f, all_con, x0, mu, epsilon, delta, Lambda, [], [], []);
% tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
% time_exact_penalty = timeit(tl1)
% time_fmincon = timeit(tmlab)

terminate_cutest_problem()


