function [lower_breakpoints, upper_breakpoints, tmin] = ...
        bounds_breakpoints(x, lb, ub, d)


    if isempty(lb)
        lower_breakpoints = inf(size(x));
    else
        lower_breakpoints = (lb - x)./d;
        l_breakpoints = ~isnan(lower_breakpoints) ...
                          & (lower_breakpoints > 0 ...
                             | (lower_breakpoints == 0 & d < 0));
        lower_breakpoints(~l_breakpoints) = inf(sum(~l_breakpoints), 1);
    end
    
    if isempty(ub)
        upper_breakpoints = inf(size(x));
    else
        upper_breakpoints = (ub - x)./d;
        u_breakpoints = ~isnan(upper_breakpoints) ...
                         & (upper_breakpoints > 0 ...
                            | (upper_breakpoints == 0 & d > 0));
        upper_breakpoints(~u_breakpoints) = inf(sum(~u_breakpoints), 1);
    end    

    if nargout >= 3
        tmin = min([lower_breakpoints; upper_breakpoints]);
    end
end
