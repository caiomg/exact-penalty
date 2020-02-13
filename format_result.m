function s = format_result(result)
% FORMAT_RESULT - 
%   
    s = '|';
    s = [s, sprintf(' %8s |', result.name)];
    if ~isempty(result.fx)
        s = [s, sprintf(' fx: % +9.3g |', result.fx)];
        s = [s, sprintf(' #f: % 6d |', result.fcount)];
        s = [s, sprintf(' err: % +9.3g |', result.error_obj)];
        s = [s, sprintf(' g lag: % +9.3g |', norm(result.lgrad))];
        s = [s, sprintf(' viol: % 9.3g |', result.nphi)];
        s = [s, sprintf(' %1d |', result.kkt)];
    else
        s = [s, sprintf(' fx:           |')];
        s = [s, sprintf(' #f:        |')];
        s = [s, sprintf(' err:           |')];
        s = [s, sprintf(' g lag:           |')];
        s = [s, sprintf(' viol:           |')];
        s = [s, sprintf(' 0 |')];
    end
end
