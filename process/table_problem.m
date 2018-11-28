problems_directory = ('/home/caio/docd/exchange/my_problems/');
for k = [4, 10, 20, 30, 40]
    problem_name  = sprintf('MADSSCHJ-%d', k);
    terminate_cutest_problem(fullfile(problems_directory, problem_name));
    prob = setup_cutest_problem(problem_name, problems_directory);
    linear_con = sum(prob.linear);
    bounds = sum(prob.bl > -1e19) + sum(prob.bu < 1e19);
    fprintf(1, '    % 4d & % 5d & % 5d & % 5d \\\\\n', prob.n, prob.m - linear_con, linear_con, bounds);
end