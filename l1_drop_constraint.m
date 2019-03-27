function [Q, R, N, ind_qr, ind_eactive] = ...
    l1_drop_constraint(cmodel, Q, R, ind_qr, ind_eactive, mu, multipliers, tol)
                                                       
% I AM ASSUMING A = QR HAS THE SAME ORDENATION AS MULTIPLIERS!
% Correct this in the future!

    if nargin < 7 || isempty(tol)
        tol = eps;
    end


    ind_j = find(multipliers < -tol | multipliers > mu + tol, 1, 'last');
    mval = min(multipliers(ind_j), mu - multipliers(ind_j));

    %[mval, ind_j] = min(min(multipliers, mu - multipliers));

    if mval < -tol
        % Remove active constraint from set
        ind_eactive(ind_eactive == ind_qr(ind_j)) = [];
        if length(ind_eactive) > length(ind_qr)
            [Q, R, ind_qr] = ...
                update_factorization(cmodel, Q, R, ind_eactive, true);
        else
            ind_qr(ind_j) = [];
            if ~isempty(ind_qr)
                [Q, R] = qrdelete(Q, R, ind_j);
                r_columns = size(R, 2);
                N = Q(:, r_columns+1:end);
            else
               R = zeros(length(Q), 0);
               Q = eye(length(Q));
               N = Q;
            end
        end
    else
        warning('cmg:runtime', 'Shouldnt have tried to drop constraint');
    end

end
