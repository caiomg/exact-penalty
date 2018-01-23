function prob = setup_cutest_problem(name, directory)
global problem_data_cutest
global problem_name_cutest
global problem_path_cutest
global all_cutest_constraints

if nargin <= 1 || isempty(directory)
    directory = '.';
end


problem_name_cutest = name;
problem_path_cutest = fullfile(directory, problem_name_cutest);

original_dir = cd(problem_path_cutest);
% Making sure we have the absolute path
problem_path_cutest = pwd();
problem_data_cutest = cutest_setup();
all_cutest_constraints = setup_cutest_constraints()
cd(original_dir);
addpath(problem_path_cutest); % I should be sanitizing this...
if nargout >= 1
    prob = problem_data_cutest;
end
end
