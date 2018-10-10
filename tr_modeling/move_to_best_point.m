function [model] = move_to_best_point(model, f)
%MOVE_TO_BEST_POINT Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 2
        f = @(x) x(1);
    end
    [dim, points_num] = size(model.points_abs);
    min_f = inf;
    min_i = 0;
    for k = 1:points_num
        val = f(model.fvalues(:, k));
        if val < min_f
            min_f = val;
            min_i = k;
        end
    end
    if min_i ~= model.tr_center
        model.tr_center = min_i;
    end
    % Here should rebuild polynomials!!!
    

end

