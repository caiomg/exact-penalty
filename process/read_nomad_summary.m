function [f_count, x] = read_nomad_summary(filename)

%     fcount_str = 'objective functions\s*=';
    fcount_str = 'blackbox evaluations                     :';
%     fval_str = 'Final f\s*=';
    x_str = 'best feasible solution                   : (';
    x_infeasible = 'best infeasible solution (min. violation): (';
    contents = fileread(filename);



    count_p = regexp(contents, fcount_str, 'end');
    [f_count, success] = sscanf(contents(count_p+1:end), '%f');
    if success < 1
        f_count = [];
        x = [];
    else
        x_pos = strfind(contents(count_p+1:end), x_str);
        if ~isempty(x_pos)
            x_pos = count_p + max(x_pos) + length(x_str);
            [x, x_read] = sscanf(contents(x_pos:end), '%f'); 
        end
        if isempty(x_pos) || x_read < 1
            x_pos = strfind(contents(count_p+1:end), x_infeasible);
            x_pos = count_p + max(x_pos) + length(x_str);
            [x, x_read] = sscanf(contents(x_pos:end), '%f');
        end
    end

end
