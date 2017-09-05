function [epsilon, Lambda, N, Q, R] = ...
    l1_criticality_step(epsilon, Lambda, current_constraints, gfx, mu, ...
                        Q, R, ind_eactive, tol_g, tol_con)

while true
    epsilon = epsilon/2;
    Lambda = Lambda/2;
    [N, Q, R, ind_eactive, ind_eviolated, ind_qr] = ...
        identify_new_constraints(current_constraints, epsilon, Q, R, ind_eactive);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_eviolated);
    if norm(N'*pseudo_gradient) > tol_g
        break;
    elseif epsilon < tol_con
        break;
    end
end

end