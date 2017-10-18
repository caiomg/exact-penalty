dim = 17;

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

all_con = cell(dim, 1);
for k = 1:dim
    Hc = zeros(dim);
    Hc(k, k) = -2;
    gc = zeros(dim, 1);
    gc(k) = 1;
    cc = (1 - round(100*(k-1)/dim)/100)^2 - 0.25;
    gk = @(x) quadratic(Hc, gc, cc, x);
    all_con{k} = gk;
end


% Initial point
x0 = -2*ones(dim, 1);
x0 = -(1+3*sqrt(2))/2*ones(dim, 1);
x0 = -1.5*ones(dim, 1);
x0 = (1:dim)';
x0 = -(1:dim)';

% Parameters
mu = 1e-2;
epsilon = 2;
delta = 1e-6;
Lambda = 0.1;

nlcon = @(x) constraints(all_con, {}, x, -1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);
x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options)
% p = @(x) l1_function(f, all_con, mu, x);
% x2_fmincon = fmincon(p, x0,[],[],[],[],[],[], nlcon, fmincon_options)

%%
x = l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda)


tl1 = @() l1_penalty(f, all_con, x0, mu, epsilon, delta, Lambda);
tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
timeit(tl1)
timeit(tmlab)




