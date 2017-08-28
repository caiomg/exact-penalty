function [h, sigma, grad_phi_j] = l1_drop_constraint(Q, R, multipliers, mu)

% I AM ASSUMING A = QR HAS THE SAME ORDENATION AS MULTIPLIERS!
% Correct this in the future!

% Hard-coded tolerance for now
tol = sqrt(eps);

if min(multipliers) < 0
    sigma = 1;
    [~, ind_j] = min(multipliers);
    grad_phi_j = Q*R(:, ind_j);
    [Q, R] = qrdelete(Q, R, ind_j);
elseif max(multipliers) > 1/mu
    sigma = -1;
    [~, ind_j] = max(multipliers);
    grad_phi_j = Q*R(:, ind_j);
    [Q, R] = qrdelete(Q, R, ind_j);
end

ind_null = sum(abs(R'), 1) < tol;
N = Q(:, ind_null);

h = sigma*(N*N')*grad_phi_j;




end