function print_results(results, fd)
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    fprintf(fd, '|');
    fprintf(fd, ' %8s |', results.name);
    if ~isempty(results.fx)
        fprintf(fd, ' f(x): % +9.3g |', results.fx);
        fprintf(fd, ' f count: % 6d |', results.fcount);
        fprintf(fd, ' f err: % +9.3g |', results.error_obj);
        fprintf(fd, ' viol.: % 9.3g |', results.nphi);
    else
        fprintf(fd, ' f(x):           |');
        fprintf(fd, ' f count:        |');
        fprintf(fd, ' f err:           |');
        fprintf(fd, ' viol.:           |');
    end
    fprintf(fd, '\n');
end
    