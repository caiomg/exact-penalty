function measure = my_measure_m_file(params)
% MY_MEASURE_M_FILE - 
%   

    p_indices = (16:27)';
    tmax = 1e4;
    computing_time = harder_problems(params, p_indices);
    S = profile_integral(computing_time, p_indices, tmax)/tmax;
    measure = 1 - S;

end

