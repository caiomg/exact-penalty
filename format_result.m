function s = format_result(result)
% FORMAT_RESULT - 
%   
    s = '|';
    s = [s, sprintf(' %8s |', result.name)];
    if ~isempty(result.fx)
        s = [s, sprintf(' f(x): % +9.3g |', result.fx)];
        s = [s, sprintf(' evals: % 6d |', result.fcount)];
        s = [s, sprintf(' f err: % +9.3g |', result.error_obj)];
        s = [s, sprintf(' g lag: % +9.3g |', norm(result.lgrad))];
        s = [s, sprintf(' viol.: % 9.3g |', result.nphi)];
        s = [s, sprintf(' %1d |', result.kkt)];
    else
        s = [s, sprintf(' f(x):           |')];
        s = [s, sprintf(' evals:        |')];
        s = [s, sprintf(' f err:           |')];
        s = [s, sprintf(' g lag:           |')];
        s = [s, sprintf(' viol.:           |')];
        s = [s, sprintf(' 0 |')];
    end
end
