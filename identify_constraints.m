function [N, Q, R, ind_eactive, ind_eviolated] = ...
                   identify_constraints(current_constraints, epsilon, ...
                   Q, R, ind_eactive)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ind_eactive_old = ind_eactive;
ind_eactive = zeros(0, 1);
ind_eviolated = zeros(0, 1);
n_constraints = size(current_constraints, 1);
for n = 1:n_constraints
    if abs(current_constraints(n).c) <= epsilon
        ind_eactive(end+1, 1) = n;
    elseif current_constraints(n).c < epsilon
        ind_eviolated(end+1, 1) = n;
    end
end

n_variables = size(Q, 1);
n_eactive = size(ind_eactive, 1);
A = zeros(n_variables, n_eactive);
for n = ind_eactive'
    A(:, n) = current_constraints(n).g;
end
% TODO: use information of previous QR decomposition
[Q, R] = qr(A);
% TODO: use QR decomposition to calculate nullspace
N = null(A');


end

