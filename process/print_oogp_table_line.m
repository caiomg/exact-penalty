function print_oogp_table_line(scenario, ipoint, mu, qo, qg, zco2, evaluations)
%PRINT_OOGP_TABLE_LINE Summary of this function goes here
%   Detailed explanation goes here

fprintf(1, '    % 2d,    % 4s,    % 10s,    % 8g,    % 8g,    % 8g    % 6d\n', mu, scenario, ipoint, qo, qg, zco2, evaluations);

end

