
                    
                    

bl = [11.570;   0.040;  17.360;   0.040;  19.600;   0.033;   5.790;  34.720];
bu = [34.720;   0.100;  38.298;   0.060;  36.400;   0.055;  57.870;  96.291];

x0_base = [16.62; 0.0502;  29.46;  0.0484;  28.0;  0.04;  31.21; 74.07];
x0_pre_gas = [19.11; 0.0502; 30.61; 0.0484; 32.2; 0.04; 26.53; 62.96];
x0_pos_gas = [14.13; 0.0502; 25.04; 0.0484; 23.8; 0.04; 35.89; 85.18];
x0_more_co2 = [16.62; 0.0577;  29.46;  0.0557;  28.0;  0.046;  31.21; 74.07];
x0_less_co2 = [16.62; 0.0427;  29.46;  0.0411;  28.0;  0.034;  31.21; 74.07];

switch initial_condition
  case 'base'
    x0 = x0_base;
  case 'pre_gas'
    x0 = x0_pre_gas;
  case 'pos_gas'
    x0 = x0_pos_gas;
  case 'more_co2'
    x0 = x0_more_co2;
  case 'less_co2'
    x0 = x0_less_co2;
  otherwise
        error('cmg:invalid_point', 'Invalid scenario choice');
end

switch scenario
  case 'C1'
    gas_lim = 200;
    co2_lim = 0.025;
  case 'C2'
    gas_lim = 130;
    co2_lim = 0.025;
  case 'C3'
    gas_lim = 200;
    co2_lim = 0.0125;
  otherwise
    error('cmg:invalid_scenario', 'Invalid scenario choice');
end



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
%

l1_options = struct('tol_radius', 1e-3, 'tol_f', 1e-6, 'eps_c', 1e-5, ...
                    'eta_1', 0, 'eta_2', 0.1, 'gamma_inc', 2, ...
                    'gamma_dec', 0.5, 'initial_radius', 0.2, ...
                    'radius_max', 2, 'criticality_mu', 50, ...
                    'criticality_beta', 10, 'criticality_omega', ...
                    0.5, 'basis', 'diagonal hessian', 'pivot_threshold', ...
                    0.1, 'poised_radius_factor', 3, 'pivot_imp', 1.1)

% Objective: minimizing negative of oil production
f = @(x) -cache.getvalue(x, 1);


all_con = {@(x) (cache.getvalue(x, 2) - gas_lim)/gas_lim;
          @(x) (cache.getvalue(x, 3) - co2_lim)/co2_lim};

mu = 1e5;
epsilon = 0.1/n_scale;
delta = 1e-6/n_scale;
Lambda = 0.1/n_scale;

[x, hs] = l1_penalty_solve(f, all_con, x0_scaled, mu, epsilon, delta, ...
                           Lambda, bl_scaled, bu_scaled, l1_options)
evaluations = ec.get_count();
  
filename =  fullfile(log_dir, [initial_condition, '_', scenario]);
save(filename);
