list_of_problems
directory = '~/docd/exchange/my_problems/';


nomad_results = {};
for k = 1:length(selected_problems)
    
    relative_filename = sprintf('%s_nomad_dec', selected_problems(k).name);
    full_summary_filename = fullfile(directory, selected_problems(k).name, relative_filename);

    [fcount, x] = read_nomad_summary(full_summary_filename);
    nomad_results{k}.name = selected_problems(k).name;
    nomad_results{k}.fcount = fcount;
    nomad_results{k}.x = x;
    
end