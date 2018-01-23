function [s, status] = simple_backtracking_line_search(f, x0, s0, min_decrease, factor, max_iter)

if nargin < 6 || isempty(max_iter)
    max_iter = 10;
end

status = false;
s = s0;
f0 = f(x0);
for k = 1:max_iter
    if f(x0 + s) < f0 - min_decrease
        status = true;
        break;
    end
    s = factor*s;
end

end