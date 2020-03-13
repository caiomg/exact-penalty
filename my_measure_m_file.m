function measure = my_measure_m_file(params)
% MY_MEASURE_M_FILE - 
%   

    p_indices = (19:30)';
    tmax = 1e4;
    measure = profile_integral(harder_problems(params, p_indices), ...
                               p_indices, tmax)/tmax;

end

