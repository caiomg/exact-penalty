function [q, q1, q2] = l1_criticality_measure(x, pseudo_gradient, Q, R, lb, ub, activities_values)

if isempty(lb)
    lb = -inf(size(x));
end
if isempty(ub)
    ub = inf(size(x));
end

dim = size(x, 1);
tol_ort = 1e-5;
tol_con = 1e-6;

[Q, R, bl_included, bu_included] = ...
    detect_and_include_active_bounds(Q, R, x, -pseudo_gradient, lb, ub, tol_con);

cols_r = size(R, 2);
N = Q(:, cols_r+1:end);

% Criticality measure
if isempty(N)
    q1 = 0;
else
    q1 = N'*pseudo_gradient;
end
% q1 = N'*pseudo_gradient;
if isempty(activities_values)
    q2 = 0;
else
    q2 = activities_values;
end
 q2 = 0; %test
q = sqrt(q1'*q1 + q2'*q2);

end

