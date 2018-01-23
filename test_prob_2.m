
% Objective
Hf = zeros(3);
gf = [0; 0; 1];
cf = 0;
f = @(x) quadratic(Hf, gf, cf, x);


all_con = cell(1, 1);

Hc = diag([10; 5; 0]);
gc = [1e-4; 1e-3; 1];
cc = 0;
gk = @(x) quadratic(Hc, gc, cc, x);
all_con{1} = gk;

Hc = zeros(3);
gc = [0; 0; -1];
cc = 0;
gk = @(x) quadratic(Hc, gc, cc, x);
all_con{2} = gk;

% Initial point
x0 = zeros(3, 1);
x0(3) = -2;
x0(2) = 0;
x0(1) = max(roots([5, 1e-4, x0(3)]));



% Parameters
mu = 100;
epsilon = 2;
delta = 1e-6;
Lambda = 0.1;

nlcon = @(x) constraints(all_con, {}, x, 1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);

%%
[x, hs2] = l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda)


tl1 = @() l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda);

time_exact_penalty = timeit(tl1)

fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);
x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options)
tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
time_fmincon = timeit(tmlab)


