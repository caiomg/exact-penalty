list_of_problems
directory = '~/docd/exchange/my_problems/';


cobyla_results = {};
for k = 1:length(selected_problems)
    
    relative_filename = sprintf('%s_cobyla_dec', selected_problems(k).name);
    full_summary_filename = fullfile(directory, selected_problems(k).name, relative_filename);

    [fcount, x] = read_cobyla_summary(full_summary_filename);
    cobyla_results{k}.name = selected_problems(k).name;
    cobyla_results{k}.fcount = fcount;
    cobyla_results{k}.x = x;
    
end