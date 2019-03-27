function [c, gc, Hc] = get_cutest_constraint(x, m, factor)

    if nargin < 3 || isempty(factor)
        factor = 1;
    end

    if nargout <= 1
        c = cutest_cons(x, m);
    else
        [c, gc] = cutest_cons(x, m);
    end
    if nargout >= 3
        Hc = cutest_ihess(x, m);
    end

    c = factor*c;
    if nargout >= 2
        gc = factor*gc;
    end
    if nargout >= 3
        Hc = factor*Hc;
    end

end
