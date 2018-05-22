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
    end
    su = sum(abs(U), 2);
    for n = 1:size(U, 1)
        if su(n) < tol
            U(n, n) = rand()*10*tol;
            degenerate = true;
        end
    end
    if degenerate
        A = U'*L';
    end


    % TODO: use information of previous QR decomposition
    [Q, R] = qr(A);
    % TODO: use QR decomposition to calculate nullspace
    N = null(A');
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
    N = null(A');
end

ind_null = sum(abs(R'), 1) < 1e-10;
N1 = Q(:, ind_null);

rank_n = size(N, 2);
if rank_n ~= size(N1, 2) || rank([N, N1], 1e-8) ~= rank_n
    1;
%     error('cmg:badnullspacerank', 'Error calculating nullspace');
else
    N = N1;
end

if size(ind_qr, 1) ~= size(R, 2)
    error('cmg:runtime_error', 'ERROR');
end


end

