function terminate_cutest_problem()

global problem_name_cutest
global problem_path_cutest
global problem_data_cutest

if ~isempty(problem_path_cutest)
    original_dir = pwd();
    cd(problem_path_cutest);
    try
        cutest_terminate();
    catch exception
        % Ignoring one exception
        if ~startsWith(exception.message, 'cutest_setup must be called first')
            retrhow(exception);
        end
    end
    cd(original_dir);
    rmpath(problem_path_cutest);

    clear global problem_path_cutest problem_name_cutest problem_data_cutest
end

end
