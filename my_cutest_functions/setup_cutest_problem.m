function prob = setup_cutest_problem(name, directory)
global problem_data_cutest
global problem_name_cutest
global problem_path_cutest
global all_cutest_constraints

if nargin <= 1 || isempty(directory)
    directory = '.';
end
if ~endsWith(directory, '/')
    % I'm using Unix style path!!!
    directory = [directory, '/'];
end

problem_name_cutest = name;
problem_path_cutest = [directory, problem_name_cutest];

original_dir = pwd();
cd(problem_path_cutest);
problem_data_cutest = cutest_setup();
all_cutest_constraints = setup_cutest_constraints()
cd(original_dir);
addpath(problem_path_cutest);
if nargout >= 1
    prob = problem_data_cutest;
end
end
