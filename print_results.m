function print_results(results, fd)
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    fprintf(fd, '|');
    fprintf(fd, ' %8s |', results.name);
    if ~isempty(results.fx)
        fprintf(fd, ' f(x): % +9.3g |', results.fx);
        fprintf(fd, ' evals: % 6d |', results.fcount);
        fprintf(fd, ' f err: % +9.3g |', results.error_obj);
        fprintf(fd, ' g lag: % +9.3g |', norm(results.lgrad));
        fprintf(fd, ' viol.: % 9.3g |', results.nphi);
        fprintf(fd, ' %1d |', results.kkt);
    else
        fprintf(fd, ' f(x):           |');
        fprintf(fd, ' evals:        |');
        fprintf(fd, ' f err:           |');
        fprintf(fd, ' g lag:           |');
        fprintf(fd, ' viol.:           |');
        fprintf(fd, ' 0 |');
    end
    fprintf(fd, '\n');
end
    