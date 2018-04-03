function [h, sigma, grad_phi_j, Q, R, ind_j] = l1_drop_constraint(Q, R, ...
                                                           multipliers, mu, tol)

% I AM ASSUMING A = QR HAS THE SAME ORDENATION AS MULTIPLIERS!
% Correct this in the future!



% TODO: try other choices
if min(multipliers) < 0
    sigma = -1;
    [~, ind_j] = min(multipliers);
    grad_phi_j = Q*R(:, ind_j);
    [Q, R] = qrdelete(Q, R, ind_j);
elseif max(multipliers) > mu
    sigma = 1;
    [~, ind_j] = max(multipliers);
    grad_phi_j = Q*R(:, ind_j);
    [Q, R] = qrdelete(Q, R, ind_j);
end

ind_null = sum(abs(R'), 1) < tol;

N1 = Q(:, ind_null);
N = null(R'*Q');
rank_n = rank(N);
if rank_n ~= rank(N1) || rank([N, N1], 1e-8) ~= rank_n
    1;
%    error('cmg:badnullspacerank', 'Error calculating nullspace');
end

h = sigma*(N*(N'*grad_phi_j));




end