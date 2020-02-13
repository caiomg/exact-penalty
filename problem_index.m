function n = problem_index(problems, problem_name)
% PROBLEM_INDEX - 
%   
    for k = 1:numel(problems)
        if strcmp(problems(k).name, problem_name)
            n = k;
            break
        end
    end
end

