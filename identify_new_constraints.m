function [N, Q, R, ind_eactive, ind_eviolated, ind_qr] = ...
                     identify_new_constraints(current_constraints, epsilon, ...
                                              Q, R, ind_eactive)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% ind_eactive = [];
ind_eviolated = zeros(0, 1);
ind_qr = zeros(0, 1);
n_constraints = size(current_constraints, 1);
for n = 1:n_constraints
    if isempty(find(ind_eactive == n, 1))
        if abs(current_constraints(n).c) <= epsilon
            ind_eactive(end+1, 1) = n;
        % Should this be here?
        % Maybe one constraint could be active and included as violated
        elseif current_constraints(n).c < epsilon
            ind_eviolated(end+1, 1) = n;
        end
    end
end

n_variables = size(Q, 1);
n_eactive = size(ind_eactive, 1);
%A = zeros(n_variables, n_eactive);
if ~isempty(ind_eactive)
    % Should check if not zero
    A = current_constraints(ind_eactive(1)).g;
    ind_qr = ind_eactive(1);
    if length(ind_eactive) > 1
        for n = ind_eactive(2:end)'
            g = current_constraints(n).g;
            if rank([A, g], 1e-8) > rank(A)
                A = [A, g];
                ind_qr = [ind_qr; n];
            end
        end
    end
    % TODO: use information of previous QR decomposition
    [Q, R] = qr(A);
    % TODO: use QR decomposition to calculate nullspace
    N = null(A');

    ind_null = sum(abs(R'), 1) < 1e-10;
    N1 = Q(:, ind_null);
    rank_n = rank(N);
    if rank_n ~= rank(N1) || rank([N, N1], 1e-8) ~= rank_n
        error('cmg:badnullspacerank', 'Error calculating nullspace');
    else
        N = N1;
    end

else
    % No active constraint
    % Test
    Q = zeros(size(Q, 1), 0);
    R = zeros(0, 0);
    N = eye(size(Q, 1));
end

end

