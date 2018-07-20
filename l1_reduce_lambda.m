function [epsilon, Lambda] = ...
    l1_reduce_lambda(epsilon, Lambda, current_constraints, gfx, mu, ...
                        Q, R, tol_g, tol_con)

n_variables = size(gfx, 1);
while true
    epsilon = epsilon/2;
    Lambda = Lambda/2;
    Q = zeros(n_variables);
    R = zeros(n_variables, 0);
    ind_eactive = zeros(0, 1);
    [ind_eactive, ind_eviolated] = ...
        identify_new_constraints(current_constraints, epsilon, ind_eactive);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                             Q, R, ind_eactive, true);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_qr, true);
    if norm(N'*pseudo_gradient) > tol_g
        break
    elseif epsilon < tol_con && Lambda < tol_g
        break
    end
end
