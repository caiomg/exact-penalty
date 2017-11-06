
terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest

% problem_name = 'WOMFLET';
problem_name = 'POLAK1';
% problem_name = 'POLAK3'; % Second derivative of cons. too much important
problem_name = 'ZY2';
prob = setup_cutest_problem(problem_name, '../my_problems/');

% Objective
f = @(x) get_cutest_objective(x);

% Constraints
n_constraints = get_cutest_total_number_of_constraints();
all_con = cell(n_constraints, 1);
for k = 1:n_constraints
    gk = @(x) get_cutest_constraint(x, k, 1);
    all_con{k} = gk;
end

% Initial point
x0 = prob.x;

% Parameters
mu = 100;
epsilon = 2;
delta = 1e-6;
Lambda = 0.1;


nlcon = @(x) constraints(all_con, {}, x, 1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);
x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options)

%%
[x, hs2] = l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda)


tl1 = @() l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda);
tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
time_exact_penalty = timeit(tl1)
time_fmincon = timeit(tmlab)

terminate_cutest_problem()


