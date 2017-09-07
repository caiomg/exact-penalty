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

% Constraint 3
Hc = -2*eye(dim);
gc = -ones(dim, 1);
cc = 3^2 - 0.5;
g3 = @(x) quadratic(Hc, gc, cc, x);

% Initial point
x0 = -2*ones(dim, 1);

% Parameters
mu = 1e-1;
epsilon = 1e-1;
delta = 1e-6;
Lambda = 0.1;

nlcon = @(x) constraints({ @(y) -g2(y), @(y) -g3(y)}, {}, x);
fmincon_options = optimoptions(@fmincon, 'Display', 'off');
x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options)
p = @(x) l1_function(f, {g2}, mu, x);
x2_fmincon = fmincon(p, x0,[],[],[],[],[],[], nlcon, fmincon_options)

%%
x = l1_penalty(f, {g2, g3}, x0, mu, epsilon, delta, Lambda)
tl1 = @() l1_penalty(f, {g2, g, g3}, x0, mu, epsilon, delta, Lambda);


tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);

% timeit(tl1)
% timeit(tmlab)




