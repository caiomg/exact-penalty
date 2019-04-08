function [Q, R, ind_qr] = update_factorization(cmodel, Q, R, ...
                                               ind_eactive, perturb)

tol = 1e-10;

n_variables = size(cmodel(1).g, 1);
n_eactive = size(ind_eactive, 1);
ind_qr = zeros(0, 1);

Q = eye(n_variables);
R = zeros(n_variables, 0);
included = false(n_eactive, 1);
n_included = 0;

for rounds = 1:n_eactive
    n_max = find(~included, 1);
    norm_max = -1;
    c_val_max = inf;
    for n = 1:n_eactive
        if ~included(n)
            this_grad = cmodel(ind_eactive(n)).g;
            norm_grad = norm(this_grad);
            c_val = cmodel(ind_eactive(n)).c;
            if norm_grad > tol
                norm_n = norm(Q(:, n_included+1:end)'*this_grad);
            else
                norm_n = norm_grad;
            end
        else
            norm_n = -1;
        end
        if norm_n > norm_max ...
                || (norm_n == norm_max && abs(c_val) < c_val_max)
            norm_max = norm_n;
            n_max = n;
            c_val_max = abs(c_val);
        end
    end
    if norm_max > tol
        [Q, R] = qrinsert(Q, R, n_included + 1, cmodel(ind_eactive(n_max)).g);
        n_included = n_included + 1;
        included(n_max) = true;
        ind_qr(n_included, 1) = ind_eactive(n_max);
    elseif perturb && n_included < n_variables
        perturbation_size = (10*tol)/(n_variables - n_included);
        grad_perturbed = cmodel(ind_eactive(n_max)).g ...
            + sum(Q(:, n_included+1:end), 2)*perturbation_size;
        [Q, R] = qrinsert(Q, R, n_included + 1, grad_perturbed);
        n_included = n_included + 1;
        included(n_max) = true;
        ind_qr(n_included, 1) = ind_eactive(n_max);
    end
end

r_columns = size(R, 2);
N = Q(:, r_columns+1:end);


if size(ind_qr, 1) ~= size(R, 2)
    error('cmg:runtime_error', 'ERROR');
end


end

