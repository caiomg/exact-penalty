function terminate_cutest_problem(problem_path)

global problem_name_cutest
global problem_path_cutest
global problem_data_cutest

if nargin < 1
    problem_path = problem_path_cutest;
end

if ~isempty(problem_path)
    original_dir = pwd();
    cd(problem_path);
    try
        cutest_terminate();
    catch exception
        % Ignoring one exception
        if ~startsWith(exception.message, 'cutest_setup must be called first')
            retrhow(exception);
        end
    end
    cd(original_dir);
    rmpath(problem_path);

    clear global problem_path_cutest problem_name_cutest problem_data_cutest
end

end
