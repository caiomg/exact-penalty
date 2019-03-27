function [prob, prob_interface] = setup_cutest_problem(name, directory)
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
    try
        problem_data_cutest = cutest_setup();
    catch setup_error
        ind = strfind(setup_error.message, ...
                      'cutest_terminate must be called first');
        if ~isempty(ind)
            cutest_terminate();
            problem_data_cutest = cutest_setup();
        else
            rethrow(setup_error);
        end
    end
    all_cutest_constraints = setup_cutest_constraints();
    cd(original_dir);
    addpath(problem_path_cutest); % I should be sanitizing this...
    if nargout >= 1
        prob = problem_data_cutest;
    end
    if nargout >= 2
        prob_interface = my_problem_interface();
        prob_interface.directory = problem_path_cutest;
        prob_interface.n_constraints = problem_data_cutest.m;
    end
end
