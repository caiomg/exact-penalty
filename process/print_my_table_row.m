function print_my_table_row(results, fd)
    
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    
    fprintf(fd, ' % 7s &', results.name);
    fprintf(fd, ' % 7d &', results.mu);
    fprintf(fd, ' % 5d &', results.fcount);
    fprintf(fd, ' % +9.3g &', results.fx);
    fprintf(fd, ' % +9.3g &', results.error_obj);
    fprintf(fd, ' % +9.3g &', norm(results.lgrad));
    fprintf(fd, ' % +9.3g \\\\\n', results.nphi);
    
