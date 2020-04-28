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

log_dir_prefix = fullfile('..', 'logs', 'oogp', datestr(now(), ...
                                                  29));
log_dir = log_dir_prefix;
suffix = 'A';
while exist(log_dir, 'dir') == 7
    log_dir = [log_dir_prefix, '-', suffix];
    suffix = char(suffix + 1);
end
success = mkdir(log_dir)
if success
    log_dir
else
    log_dir = '.';
end

all_mu = [1024, 1000, 512, 2000];
all_scenarios = {'C3', 'C2', 'C1'};
all_ic = {'base', 'pre_gas', 'pos_gas', 'more_co2', 'less_co2'};
for mu = all_mu
    for ic = 1:length(all_ic)
        for scen = 1:length(all_scenarios)
            scenario = all_scenarios{scen}
            initial_condition = all_ic{ic}
            experiment_oogp
        end
    end
end

