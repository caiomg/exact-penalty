function [q, q1, q2] = l1_criticality_measure(x, pseudo_gradient, N, bl, bu, activities_values)

% Criticality measure
q1 = project_to_bounds(x - N*(N'*pseudo_gradient), bl, bu) - x;
% q1 = N'*pseudo_gradient;
if isempty(activities_values)
    q2 = 0;
else
    q2 = activities_values;
end
q = sqrt(q1'*q1 + q2'*q2);

end