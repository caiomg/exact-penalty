function [Q, R, N, ind_eactive] = ...
    l1_drop_constraint(Q, R, N, ind_eactive, mu, multipliers, tol)
                                                       
% I AM ASSUMING A = QR HAS THE SAME ORDENATION AS MULTIPLIERS!
% Correct this in the future!

    if nargin < 7 || isempty(tol)
        tol = eps;
    end

    rows_qr = size(Q, 2) - size(N, 2);

    % TODO: try other choices
    if min(multipliers) < tol && ~(-min(multipliers) < max(multipliers) - mu)
        [~, ind_j] = min(multipliers);
        [Q, R] = qrdelete(Q, R, ind_j);
    elseif max(multipliers) > mu + tol
        [~, ind_j] = max(multipliers);
        [Q, R] = qrdelete(Q, R, ind_j);
    end

    % Remove active constraint from set
    ind_eactive(ind_j) = [];

    % Update space orthogonal to active constraints
    rows_qr = rows_qr - 1;
    ind_null = (rows_qr+1):size(Q, 2);
    N = Q(:, ind_null);

end