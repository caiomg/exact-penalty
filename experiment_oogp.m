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
if exist('santos_gas_network', 'file') ~= 2
    addpath('../network_problem/')
end

l1_options = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_1', 0, 'eta_2', 0.1, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 0.5, 'radius_max', 1e3, ...
                        'criticality_mu', 50, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'diagonal hessian', ...
                        'pivot_threshold', 0.1, 'poised_radius_factor', 2, ...
                        'pivot_imp', 1.1)
                    
                    

bl = [11.570;   0.040;  17.360;   0.040;  19.600;   0.033;   5.790;  34.720];
bu = [34.720;   0.100;  38.298;   0.060;  36.400;   0.055;  57.870;  96.291];
x0 = [16.62; 0.0426;  29.46;  0.0411;  28.0;  0.034;  31.21;  74.07];

O = (bl + bu)/2;
fixed_scale = (bu - bl)/2;
no_scale = isinf(fixed_scale);
fixed_scale(no_scale) = no_scale(no_scale);
n_scale = sqrt(norm(fixed_scale));

x0_scaled = (x0 - O)./fixed_scale;
bl_scaled = (bl - O)./fixed_scale;
bu_scaled = (bu - O)./fixed_scale;

ec = evaluation_counter(@(y) santos_gas_network(y.*fixed_scale + O));
cache = simple_cache(@(x) ec.evaluate(x), 3);

f = @(x) cache.getvalue(x, 1);
all_con = {@(x) (cache.getvalue(x, 2) - 200)/1e4;
          @(x) (cache.getvalue(x, 3) - 0.025)};

mu = 10;
epsilon = 1;
delta = 1e-6;
Lambda = 0.1;

[x, hs] = l1_penalty_solve(f, all_con, x0_scaled, mu, epsilon, delta, Lambda, bl_scaled, bu_scaled, l1_options)

  

