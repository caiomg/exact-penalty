
p_indices = (19:30)';
tmax = 1e4;
measure = @(params) profile_integral(harder_problems(params, p_indices), ...
                                     p_indices, tmax)/tmax;
params_initial = [1;
                  0.8;
                  0.01;
                  0.001];

tuned_parameters = bfo(measure, params_initial, ...
                       'xlower', [0.01; 0.001; 0.01; 0.001], ...
                       'xupper', [1000, 10, 0.5, 0.15], ...
                       'xscale', [1, 1, 0.1, 0.1])
