function print_times(x, computing_times, m)

logfile = 'x_computing_times.txt';
try
    fid = fopen(logfile, 'a');
    fprintf(fid, '%s', format_x_times(x, computing_times, m));
    fclose(fid);
catch thiserror
    'Ignore';
end


end

function s = format_x_times(x, times, m)
    s = [newline, 'X: '];
    s = [s, sprintf(' %g', x)];
    s = [s, newline];
    s = [s, 'T: '];
    s = [s, sprintf(' % 5d', times)];
    s = [s, newline];
    s = [s, sprintf('M = % 10g', m)];
    s = [s, newline];
end

    