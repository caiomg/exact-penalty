function [fx, dfx, d2fx] = quadratic(H, g, c, x)

    fx = c + g'*x + 0.5*(x'*H*x);
    dfx = g + H*x;
    d2fx = H;

end