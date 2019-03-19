function [f_count, x] = read_cobyla_summary(filename)

    magic_str = '\** CUTEst statistics \**';
%     fcount_str = 'objective functions\s*=';
    fcount_str = 'FEVALS =';
%     fval_str = 'Final f\s*=';
    x_str = 'X =';
    contents = fileread(filename);

    final_stats_idx = regexp(contents, magic_str);
    contents = contents(final_stats_idx:end);

    count_p = regexp(contents, fcount_str, 'end');
    [f_count, charcount] = sscanf(contents(count_p+1:end), '%f');

    x_pos = strfind(contents(charcount:end), x_str);
    x_pos = max(x_pos) + length(x_str);
    x = sscanf(contents(x_pos:end), '%f');

end
