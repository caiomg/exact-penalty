function [f_count, fval] = read_cutest_summary(filename)

    magic_str = '\** CUTEst statistics \**';
    fcount_str = 'objective functions\s*=';
    fval_str = 'Final f\s*=';


    contents = fileread(filename);
    final_stats_idx = regexp(contents, magic_str);

    contents = contents(final_stats_idx:end);

    count_p = regexp(contents, fcount_str, 'end');
    f_count = sscanf(contents(count_p+1:end), '%f');

    fval_p = regexp(contents, fval_str, 'end');
    fval = sscanf(contents(fval_p+1:end), '%f');

end
