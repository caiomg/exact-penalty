
terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest

problem_name = 'HS100';
% problem_name = 'WOMFLET';
% problem_name = 'POLAK1';
% problem_name = 'SNAKE';
% problem_name = 'POLAK3'; % Second derivative of cons. too much important
% % problem_name = 'QC';
problem_name = 'CB2';
% problem_name = 'LOOTSMA';
% % problem_name = 'HS88';
% problem_name = 'HS18';
% problem_name = 'HS19';
% problem_name = 'HS21';
% problem_name = 'HS101';
problem_name = 'HS101';

prob = setup_cutest_problem(problem_name, '../my_problems/');

% Objective
f_obj = @(x) get_cutest_objective(x);
counter = evaluation_counter(f_obj);
f = @(x) counter.evaluate(x);

% Constraints
n_constraints = sum(prob.cl > -1e19) + sum(prob.cu < 1e19);
dim = prob.n;

% Bound constraints
bl = [];
bu = [];
bl = prob.bl;
bu = prob.bu;
cutest_lower_bounds = prob.bl > -1e19;
cutest_upper_bounds = prob.bu < 1e19;
lower_bounds = bl > -1e19;
upper_bounds = bu < 1e19;

% 'Nonlinear' constraints
all_con = cell(prob.m, 1);

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

% Initial point
x0 = prob.x;

% Parameters
mu = 10;

epsilon = 1;
delta = 1e-6;
Lambda = 0.1;


nlcon = @(x) constraints(all_con, {}, x, 1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', false);

global x_fmincon
[x_fmincon, fx_fmincon, exitflag, output, lambda_fmincon] = fmincon(f, x0,[],[],[],[], bl, bu, nlcon, fmincon_options)
nlcon_fmincon = max(0, nlcon(x_fmincon));

fmincon_count = counter.get_count()
counter.reset_count()
% counter.set_max_count(15000);



l1_options = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_1', 0, 'eta_2', 0.05, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 1, 'radius_max', 1e3, ...
                        'criticality_mu', 50, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'diagonal hessian', ...
                        'pivot_threshold', 0.001, 'poised_radius_factor', 6, ...
                        'pivot_imp', 1.1, 'debug', true, 'inspect_iteration', 30)

%%
warning('off', 'cmg:badly_conditioned_system');
p_seed = rng('default');
% [x, hs2] = l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda)
[x, hs2] = l1_penalty_solve(f, all_con, x0, mu, epsilon, delta, Lambda, bl, bu, l1_options)
fx = f(x)
nphi = norm(max(0, nlcon(x)))
warning('on', 'cmg:badly_conditioned_system');


rng(p_seed);

l1_count = counter.get_count()
counter.reset_count()


tl1 = @() l1_penalty_solve(f, all_con, x0, mu, epsilon, delta, Lambda, [], [], []);
% tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
% time_exact_penalty = timeit(tl1)
% time_fmincon = timeit(tmlab)

[kkt_ok, lgrad] = check_kkt(f, all_con, x, bl, bu, 1e-5, 5e-5)

terminate_cutest_problem()


