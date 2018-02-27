function [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive)

tol = 1e-8;

n_variables = size(Q, 1);
n_eactive = size(ind_eactive, 1);
ind_qr = zeros(0, 1);
A = zeros(n_variables, 0);
rank_a = 0;
for n = 1:n_eactive
    if norm(current_constraints(ind_eactive(n)).g, 1) > tol
        if rank([A, current_constraints(ind_eactive(n)).g], tol) > rank_a
            A = [A, current_constraints(ind_eactive(n)).g];
            rank_a = rank_a + 1;
            ind_qr = [ind_qr; ind_eactive(n)];
        end
    end
end
% TODO: use information of previous QR decomposition
[Q, R] = qr(A);

ind_null = sum(abs(R'), 1) < 1e-10;
N1 = Q(:, ind_null);
if size(N1, 2) < size(A, 1) - rank(A)
    N = null(A');
    %error('cmg:badnullspacerank', 'Error calculating nullspace');
else
    N = N1;
end

if size(ind_qr, 1) ~= size(R, 2)
    error('cmg:runtime_error', 'ERROR');
end


end

