function [f_count, fval, x] = read_cutest_summary(filename)

    magic_str = '\** CUTEst statistics \**';
    fcount_str = 'objective functions\s*=';
    fval_str = 'Final f\s*=';
    x_str = 'X =';

    contents = fileread(filename);
    
    
    x_pos = strfind(contents, x_str);
    x_pos = max(x_pos) + length(x_str);
    x = sscanf(contents(x_pos:end), '%f');
    
    final_stats_idx = regexp(contents, magic_str);

    contents = contents(final_stats_idx:end);

    count_p = regexp(contents, fcount_str, 'end');
    f_count = sscanf(contents(count_p+1:end), '%f');

    fval_p = regexp(contents, fval_str, 'end');
    fval = sscanf(contents(fval_p+1:end), '%f');

end
