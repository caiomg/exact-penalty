function print_oogp_table_line_tex(scenario, ipoint, mu, qo, qg, zco2, evaluations)
%PRINT_OOGP_TABLE_LINE Summary of this function goes here
%   Detailed explanation goes here

fprintf(1, '    % 4s  &  % 10s  &  % 8.0f  &  % 3.1f  &  % 1.2f  &  % 6f\n', scenario, ipoint, qo/10, qg, zco2*100, evaluations);

end

