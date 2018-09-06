function [f, g, H] = scale_function(fun, O, S, w)

    if nargout <= 1
        f = fun(O + S*w);
    else
        if nargout == 3
            [f, g, H] = fun(O + S*w);
            H = S'*H*S;
        else
            [f, g] = fun(O + S*w);
        end
        g = S'*g;
    end

end
