
terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest

problem_name = 'HS100';
% problem_name = 'WOMFLET';
problem_name = 'POLAK1';
% problem_name = 'SNAKE';
problem_name = 'POLAK3'; % Second derivative of cons. too much important
% problem_name = 'QC';
% problem_name = 'CB2';
% problem_name = 'LOOTSMA';
% problem_name = 'HS88';
prob = setup_cutest_problem(problem_name, '../my_problems/');

% Objective
f_obj = @(x) get_cutest_objective(x);
counter = evaluation_counter(f_obj);
f = @(x) counter.evaluate(x);

% Constraints
n_constraints = get_cutest_total_number_of_constraints();
all_con = cell(n_constraints, 1);
for k = 1:n_constraints
    gk = @(x) evaluate_my_cutest_constraint(x, k, 1);
    all_con{k} = gk;
end

% Initial point
x0 = prob.x;

% Parameters
mu = 100;
if strcmp(problem_name, 'SNAKE')
    mu = 10000;
end
epsilon = 0.85;
delta = 1e-6;
Lambda = 0.075;


nlcon = @(x) constraints(all_con, {}, x, 1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);

x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
fx_fmincon = f(x_fmincon);

counter.get_count()
counter.reset_count()

%%
p_seed = rng('default');
% [x, hs2] = l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda)
[x, hs2] = l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda)
fx = f(x)
nphi = norm(max(0, nlcon(x)))

rng(p_seed);

counter.get_count()
counter.reset_count()


tl1 = @() l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda);
% tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
% time_exact_penalty = timeit(tl1)
% time_fmincon = timeit(tmlab)

terminate_cutest_problem()


