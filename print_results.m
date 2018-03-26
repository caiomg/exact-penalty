function print_results(results, fd)
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    fprintf(fd, '|');
    fprintf(fd, '%8s  |', results.name);
    fprintf(fd, 'f count: % 4d  |', results.fcount);
    fprintf(fd, 'f(x): % +9.3g  |', results.fx);
    fprintf(fd, 'f count ml: % 4d  |', results.fcount_fmincon);
    fprintf(fd, 'f error: % +9.3g  |', results.error_obj);
    fprintf(fd, 'x error: % +9.3g  |', results.error_x);
    fprintf(fd, '\n');
end
    