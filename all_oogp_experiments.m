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

log_dir = fullfile('..', 'logs', 'oogp', datestr(now(), 29));
all_scenarios = {'C1', 'C2' 'C3'};
all_ic = {'base', 'pre_gas', 'pos_gas', 'more_co2', 'less_co2'};
for scen = 1:length(all_scenarios)
    for ic = 1:length(all_ic)
        scenario = all_scenarios{scen}
        initial_condition = all_ic{ic}
        experiment_oogp
    end
end

