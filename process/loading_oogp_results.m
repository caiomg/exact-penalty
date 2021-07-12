
new_results = struct('qo', {}, 'qg', {}, 'zco2', {}, ...
                     'evaluations', {}, 'initial_condition', {}, 'mu', {})
scenarios = struct('name', {'C1', 'C2', 'C3'}, 'results', {new_results, new_results, new_results})

directory = '/home/caio/pdoc/exact-penalty-paper/logs/2020-04-30';
query = fullfile(directory, '*_1000.mat');

files = dir(query);


for m = 1:numel(files)
    filename = fullfile(files(m).folder, files(m).name);
    for k = 1:numel(scenarios)
        if contains(filename, scenarios(k).name)
            current_results = load(filename, ...
                'qo', 'qg', 'zco2', 'evaluations', 'initial_condition', 'mu');
            scenarios(k).results(end+1) = current_results;
            break
        end
    end
end


qo1 = [scenarios(1).results.qo]'
best1 = max(qo1)
mean1 = mean(qo1)
evals1 = [scenarios(1).results.evaluations]'

qo2 = [scenarios(2).results.qo]'
best2 = max(qo2)
mean2 = mean(qo2)
qg2 = [scenarios(2).results.qg]'
viol2 = max(0, max(qg2) - 130)
evals2 = [scenarios(2).results.evaluations]'

qo3 = [scenarios(3).results.qo]'
best3 = max(qo3)
mean3 = mean(qo3)
zco23 = [scenarios(3).results.zco2]'
viol3 = max(0, max(zco23) - 0.0125)
evals3 = [scenarios(3).results.evaluations]'
