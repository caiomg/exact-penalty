function [epsilon, Lambda, A, N, ind_eactive, ind_eviolated, pseudo_gradient] = ...
    l1_criticality_step(epsilon, Lambda, current_constraints, gfx, mu, Q, R, tol_g, tol_con)

while true
    epsilon = epsilon/2;
    Lambda = Lambda/2;
    [A, N, ind_eactive, ind_eviolated] = ...
        identify_constraints(current_constraints, epsilon, Q, R);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_eactive);
    if norm(N'*pseudo_gradient) > tol_g
        break;
    elseif epsilon < tol_con
        break;
    end
end

end