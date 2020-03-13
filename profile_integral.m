function S = profile_integral(times, indices, tmax)

    other_times = get_other_times(indices);
    computing_time = [other_times];
    [rho, tau] = dm_performance_profile(computing_time);
    [rho, tau] = adjust_performance_profile_scale(rho, tau, tmax);
    n_steps = size(rho, 1);
    S = 0;
    for k = 1:n_steps-1
        S = S + rho(k)*(tau(k+1) - tau(k));
    end
    S = S;
end


function times = get_other_times(indices)

al_evals = [805, 288, 302; % HS18
            nan, 2808, 1535; % HS19
            57, 60, 177; % HS21
            613, 504, 841; %HS23
            95, 81, 338; %HS30
            506, 277, 427; %HS31
            8505, 850, 1640; %HS34
            235, 181, 502354; %HS36
            2790, 888, 1422; %HS37
            223, 114, 196; %HS59
            9962, 976, 512; %HS65
            11783, 5609, 831; %HS66
            4473, 848, 965; %HS67
            5424, 656, 455; %HS70
            8061, 26804, 6720; %HS72
            1777, 1311, 7752; %HS83
            678500, 607443, 4152499; %HS84
            427470, nan, 60009; %HS85
            215, 377, 1936; %HS95
            223, 377, 1940; %HS96
            1754, 700, 7234029; %HS97
            1754, 700, 7234029; %HS98
            215385, 202679, 144597; %HS101
            23618, 267740, 1068121; %HS102
            54425, 222196, 1820078; %HS103
            87385, 19614, 7993; %HS104
            6891, 6193, 12149; %HS105
            657806, 37640, 504366; %HS106
            136210, 99678, 6215728; %HS116
            220633, 3633, 186851 %HS118
            ];

    times = al_evals(indices, :);
end

function [rho, tau] = adjust_performance_profile_scale(rho, tau, tmax)

    n_elements = size(rho, 1);
    assert(n_elements == numel(tau));
    
    cut_point = find(tau >= tmax, 1);
    if ~isempty(cut_point)
        if tau(cut_point) == tmax
            tau = tau(1:cut_point, :);
            rho = rho(1:cut_point, :);
        else
            tau = [tau(1:cut_point-1, :);
                   tmax];
            rho = [rho(1:cut_point-1, :);
                   rho(cut_point-1, :)];
        end
    else
        tau = [tau;
               tmax];
        rho = [rho;
               rho(end, :)];
    end
    
end
