function print_my_table_row(results, fd)
    
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    
    fprintf(fd, ' % 7s &', results.name);
    if log10(results.mu) - int64(log10(results.mu)) <= 100*eps
        fprintf(fd, ' $10^{% 2d}$ &', log10(results.mu));
    else
        fprintf(fd, ' % 8d &', results.mu);
    end
    fprintf(fd, ' % 6d &', results.fcount);
    if isfield(results, 'best_fcount')
        fprintf(fd, ' % 6d &', results.best_fcount);
    end
    fprintf(fd, ' % +9.3g &', results.fx);
    fprintf(fd, ' % 9.3g &', norm(results.lgrad));
    fprintf(fd, ' % 9.3g \\\\\n', results.nphi);
    
