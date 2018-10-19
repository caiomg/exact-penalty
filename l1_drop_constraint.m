function [Q, R, N, ind_qr, ind_eactive] = ...
    l1_drop_constraint(cmodel, Q, R, ind_qr, ind_eactive, mu, multipliers, tol)
                                                       
% I AM ASSUMING A = QR HAS THE SAME ORDENATION AS MULTIPLIERS!
% Correct this in the future!

    if nargin < 7 || isempty(tol)
        tol = eps;
    end


    % TODO: try other choices
    if min(multipliers) < tol && ~(-min(multipliers) < max(multipliers) - mu)
        [~, ind_j] = min(multipliers);
    elseif max(multipliers) > mu + tol
        [~, ind_j] = max(multipliers);
    end

    % Remove active constraint from set
    ind_eactive(ind_eactive == ind_qr(ind_j)) = [];
    if length(ind_eactive) > length(ind_qr)
        [N, Q, R, ind_qr] = update_factorization(cmodel, ...
                                                 Q, R, ind_eactive, true);
    else
        ind_qr(ind_j) = [];
        if length(ind_qr) > 0
            [Q, R] = qrdelete(Q, R, ind_j);
            % Update space orthogonal to active constraints
            r_columns = size(R, 2);
            N = Q(:, r_columns+1:end);
        else
           R = zeros(length(Q), 0);
           Q = eye(length(Q));
           N = Q;
        end
    end

end