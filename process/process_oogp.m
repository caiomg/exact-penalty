if exist('santos_gas_network', 'file') ~= 2
    addpath('../network_problem/')
end



all_mu = [6, 3, 4, 5];
all_scenarios = {'C3', 'C2', 'C1'};
all_ic = {'base', 'pre_gas', 'pos_gas', 'more_co2', 'less_co2'};
log_dir = '/home/caio/docd/exchange/logs/oogp/2018-11-04-B';
for muexp = all_mu
    for scen = 1:length(all_scenarios)
    scenario = all_scenarios{scen};
        for ic = 1:length(all_ic)
            initial_condition = all_ic{ic};
            rel_name = sprintf('%s_%s_%d.mat', initial_condition, scenario, muexp);
            filename =  fullfile(log_dir, rel_name);
            if exist(filename, 'file') == 2
                V = load(filename, 'x_full_space', 'evaluations');
                [qo, qg, zco2] = santos_gas_network(V.x_full_space);
                print_oogp_table_line(scenario, initial_condition, muexp, qo, qg, zco2, V.evaluations);
            end
        end
    end
end
