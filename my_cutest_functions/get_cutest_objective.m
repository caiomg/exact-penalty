function [f, g, H] = get_cutest_objective(x, factor)

    if nargin < 2 || isempty(factor)
        factor = 1;
    end

    if nargout == 1
        f = cutest_obj(x);
    else
        [f, g] = cutest_obj(x);
    end
    if nargout >= 3
        H = cutest_ihess(x, 0);
    end

    f = factor*f;
    if nargout >= 2
        g = factor*g;
    end
    if nargout >= 3
        H = factor*H;
    end

end
