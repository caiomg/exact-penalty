function [epsilon, Lambda] = ...
    l1_reduce_lambda(epsilon, Lambda, current_constraints, gfx, mu, ...
                        Q, R, x, bl, bu, tol_g, tol_con)

n_variables = size(gfx, 1);
while true
%     epsilon = epsilon/2;
    Lambda = Lambda/2;
    Q = zeros(n_variables);
    R = zeros(n_variables, 0);
    ind_eactive = zeros(0, 1);
    [ind_eactive, ind_eviolated] = ...
        identify_new_constraints(current_constraints, epsilon, ind_eactive);
    [Q, R, ind_qr] = update_factorization(current_constraints, ...
                                             Q, R, ind_eactive, true);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_qr, true);
    measure = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, ...
                                    bu, []);
    if measure > tol_g
        break
    elseif Lambda < tol_g
        break
    end
end
