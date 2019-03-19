function t = tr_radius_breakpoint(d, radius, s0)

    if nargin < 3
        t = radius/norm(d);
    else
        t = roots([d'*d, 2*s0'*d, s0'*s0 - radius^2]);
        t = min(t(t > 0));
        if isempty(t) || ~isreal(t)
            t = radius/norm(d);
        end
    end
end