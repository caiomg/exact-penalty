function [Q, R, ind_qr, ind_eactive, ind_dropped, multiplier_dropped] = ...
    l1_drop_constraint(cmodel, Q, R, ind_qr, ind_eactive, mu, multipliers, tol)
                                                       
% I AM ASSUMING A = QR HAS THE SAME ORDENATION AS MULTIPLIERS!
% Correct this in the future!

    if nargin < 7 || isempty(tol)
        tol = eps;
    end


%     ind_j = find(multipliers < -tol | multipliers > mu + tol, 1, 'last');
%     mval = min(multipliers(ind_j), mu - multipliers(ind_j));

    %[mval, ind_j] = min(min(multipliers, mu - multipliers));

    c_vals = [cmodel(ind_qr).c]';

    best_candidates = (multipliers < -tol & c_vals <= 0) ...
                       | (multipliers > mu + tol & c_vals >= 0);
    chosen_c_absval = -inf;
    for m = 1:length(ind_qr)
       if best_candidates(m) ...
               && chosen_c_absval < abs(c_vals(m))
           ind_j = m;
           chosen_c_absval = abs(c_vals(m));
       end
    end
    if isinf(chosen_c_absval)
        for m = 1:length(ind_qr)
            if (multipliers(m) < -tol ...
                || multipliers(m) > mu + tol) ...
                    && chosen_c_absval < abs(c_vals(m))
                ind_j = m;
                chosen_c_absval = abs(c_vals(m));
            end
        end
    end
    
    if ~exist('ind_j', 'var')
        warning();
    end

%    ind_j = find(multipliers < -tol | multipliers > mu + tol, 1, 'last');
    mval = min(multipliers(ind_j), mu - multipliers(ind_j));



    if mval < -tol
        % Remove active constraint from set
        ind_dropped = ind_qr(ind_j);
        multiplier_dropped = multipliers(ind_j);
        ind_eactive(ind_eactive == ind_dropped) = [];
        if length(ind_eactive) >= length(ind_qr)
            [Q, R, ind_qr] = ...
                update_factorization(cmodel, Q, R, ind_eactive, true);
        else
            ind_qr(ind_j) = [];
            if ~isempty(ind_qr)
                [Q, R] = qrdelete(Q, R, ind_j);
            else
               R = zeros(length(Q), 0);
               Q = eye(length(Q));
            end
        end
    else
        ind_dropped = [];
        multiplier_dropped = [];
        warning('cmg:runtime', 'Shouldnt have tried to drop constraint');
    end

end
