function [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive, perturb)

tol = 1e-8;

n_variables = size(current_constraints(1).g, 1);
n_eactive = size(ind_eactive, 1);
ind_qr = zeros(0, 1);
A = zeros(n_variables, 0);
if perturb
    for n = 1:n_eactive
        norm_n = norm(current_constraints(ind_eactive(n)).g, 1);
        if norm_n > tol
            A = [A, current_constraints(ind_eactive(n)).g];
            ind_qr(end+1, 1) = ind_eactive(n);
        end
    end
    % Testing degeneracy
    [L, U, P] = lu(A');
    degenerate = false;
    if size(L, 1) > n_variables
        L = L(1:n_variables, :);
        ind_qr = P*ind_qr;
        ind_qr = ind_qr(1:n_variables, :);
        degenerate = true;
        A = U'*L';
    end
    [Q, R] = qr(A);
    for n = 2:size(R, 2)
        if abs(R(n, n)) < tol
            R(n, n) = max(1, max(diag(R)))*tol*10*rand();
            degenerate = true;
        end
    end
    A = Q*R;

else
    Q = eye(n_variables);
    R = zeros(n_variables, 0);
    inserted = 0;
    for n = 1:n_eactive
        this_grad = current_constraints(ind_eactive(n)).g;
        norm_n = norm(A*(R\(Q'*this_grad)) - this_grad, 1);
        if norm_n > tol
            % Try to add column
            [Q, R] = qrinsert(Q, R, inserted + 1, this_grad);
            % Check if column really added
            if size(R, 2) > inserted
                inserted = inserted + 1;
                A = [A, this_grad];
                ind_qr(end+1, 1) = ind_eactive(n);
            end
        end
    end
end

r_columns = size(R, 2);
N = Q(:, r_columns+1:end);


if size(ind_qr, 1) ~= size(R, 2)
    error('cmg:runtime_error', 'ERROR');
end


end

