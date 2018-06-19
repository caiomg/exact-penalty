function print_results(results, fd)
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    fprintf(fd, '|');
    fprintf(fd, '%8s  |', results.name);
    fprintf(fd, 'f(x): % +9.3g  |', results.fx);
    fprintf(fd, 'f count: % 6d  |', results.fcount);
    fprintf(fd, 'f error: % +9.3g  |', results.error_obj);
    fprintf(fd, 'violation: % +9.3g  |', results.nphi);
    fprintf(fd, '\n');
end
    