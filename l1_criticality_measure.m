function [q, q1, q2] = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, activities_values)

if isempty(bl)
    bl = -inf(size(x));
end
if isempty(bu)
    bu = inf(size(x));
end

dimension = size(x, 1);
tol_ort = 1e-5;
tol_con = 1e-10;
while true
    cols_r = size(R, 2);
    N = Q(:, cols_r+1:end);
	% Test bounds
    l_active = x - bl <= tol_con & N*(N'*pseudo_gradient) > 0;
    u_active = x - bu >= -tol_con & N*(N'*pseudo_gradient) < 0;
    n_lower_active = sum(l_active);
    n_upper_active = sum(u_active);
    I = eye(dimension);
    Al = I(:, l_active);
    Au = -I(:, u_active);
    bound_included = false;
    for k = 1:n_lower_active
        this_grad = Al(:, k);
        norm_n = norm((Q*R)*(R\(Q'*this_grad)) - this_grad, 1);
        if norm_n > tol_ort
            [Q, R] = qrinsert(Q, R, 1, this_grad);
            bound_included = true;
            break
        else
            l_active(this_grad ~= 0) = false;
        end
    end
    if bound_included
        continue
    end
    for k = 1:n_upper_active
        this_grad = Au(:, k);
        norm_n = norm((Q*R)*(R\(Q'*this_grad)) - this_grad, 1);
        if norm_n > tol_ort
            [Q, R] = qrinsert(Q, R, 1, this_grad);
            bound_included = true;
            break
        else
            u_active(this_grad ~= 0) = false;
        end
    end
    if ~bound_included
        break
    end
end

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

