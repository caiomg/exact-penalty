function [f_count, fval] = read_cutest_summary(filename)

    magic_str = ['************************ CUTEst statistics **********' ...
                 '**************'];
    fcount_str = '# objective functions   =';
    fval_str = 'Final f                 =';


    contents = fileread(filename);
    final_stats_idx = strfind(contents, magic_str);

    contents = contents(final_stats_idx:end);
    count_p = strfind(contents, fcount_str);
    contents = contents(count_p + length(fcount_str):end);
    f_count = sscanf(contents, '%f');

    fval_p = strfind(contents, fval_str);
    contents = contents(fval_p + length(fval_str):end);
    fval = sscanf(contents, '%f');

end
