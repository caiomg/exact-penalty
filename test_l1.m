dim = 2;

% Objective
Hf = 2*eye(dim);
gf = -4*ones(dim, 1);
cf = dim*4;
f = @(x) quadratic(Hf, gf, cf, x);

% Constraint
Hc = -2*eye(dim);
gc = zeros(dim, 1);
cc = 1;
g = @(x) quadratic(Hc, gc, cc, x);

% Constraint 2
Hc = -2*eye(dim);
gc = zeros(dim, 1);
gc(1) = 2;
cc = (0.5)^2 - 1;
g2 = @(x) quadratic(Hc, gc, cc, x);

% Initial point
x0 = -2*ones(dim, 1);

% Parameters
mu = 0.2;
epsilon = 1;
delta = 0.1;
Lambda = 1;

x = l1_penalty(f, {g2, g}, x0, mu, epsilon, delta, Lambda)
tl1 = @() l1_penalty(f, {g2, g}, x0, mu, epsilon, delta, Lambda);

nlcon = @(x) constraints({@(y) -g(y), @(y) -g2(y)}, {}, x);
fmincon_options = optimoptions(@fmincon, 'Display', 'off');
x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options)
tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);

timeit(tl1)
timeit(tmlab)