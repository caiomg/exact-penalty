function x = l1_penalty(f, phi, x0, mu, epsilon, delta, Lambda )
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here


x = x0;
n_variables = size(x, 1);
n_constraints = size(phi, 1);
tol_g = n_variables*sqrt(eps);
tol_con = n_constraints*sqrt(eps);

p = @(x) l1_function(f, phi, mu, x);
% QR decomposition of constraints gradients matrix A
Q = zeros(n_variables, 0);
R = zeros(0, 0);
ind_eactive = zeros(0, 1);



current_constraints = evaluate_constraints(phi, x);

while true
    [N, Q, R, ind_eactive, ind_eviolated] = ...
                      identify_constraints(current_constraints, epsilon, ...
                                           Q, R, ind_eactive);

    [fx, gfx, Hfx] = f(x);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_eviolated);

    if (norm(N'*pseudo_gradient) > Lambda)
        B = zeros(size(Hfx));
        for n = ind_eviolated'
            B = B - (current_constraints(n).H)/mu;
        end
        B = B + Hfx;
        % Calculate Newton direction
        u = (N'*B*N)\(N'*pseudo_gradient);
        h = -N*u;
        x = l1_linear_search(pseudo_gradient, h, x, current_constraints, ...
                             p, epsilon, delta);
        current_constraints = evaluate_constraints(phi, x);
    elseif (~isempty(ind_eviolated) || (epsilon > tol_con) || ...
            (norm(N'*pseudo_gradient) > tol_g))
        % calculate multipliers
        multipliers = R\(Q'*pseudo_gradient);

        % Are there conditions for dropping one constraint?
        if sum(multipliers < 0 | 1/mu < multipliers)
            [h, sigma, grad_phi_j] = l1_drop_constraint(Q, R, multipliers, mu);
            if ((pseudo_gradient + min(0, sigma)*grad_phi_j)'*h < -delta)
                x = l1_linear_search(pseudo_gradient, h, x, ...
                                     current_constraints, p, epsilon, delta);
                current_constraints = evaluate_constraints(phi, x);
            else
                [epsilon, Lambda] = l1_criticality_step(epsilon, Lambda, ...
                                                        current_constraints,...
                                                        gfx, mu, Q, R, ...
                                                        ind_eactive, ...
                                                        tol_g, tol_con);
            end
        else
            B = zeros(size(Hfx));
            for n = ind_eviolated'
                B = B - (current_constraints(n).H)/mu;
            end
            for n = ind_eactive'
                B = B - multipliers(n)*(current_constraints(n).H);
            end
            B = B + Hfx;
            % Calculate Newton direction
            u = (N'*B*N)\(N'*pseudo_gradient);
            h = -N*u;
            % Recalculate constraints
            phih = zeros(size(ind_eactive));
            for n = 1:ind_eactive'
               phih(n) = phi{n}(x + h);
            end
            v = Q*(R'\phih);
            normphi = norm([current_constraints(ind_eactive).c]);
            ppgrad = N'*pseudo_gradient;
            if (p(x + h + v) <= p(x) - delta*(norm(ppgrad)^2 + normphi))
                x = x + h + v;
                current_constraints = evaluate_constraints(phi, x);
            else
                [epsilon, Lambda] = l1_criticality_step(epsilon, Lambda, ...
                                                        current_constraints, ...
                                                        gfx, mu, Q, R, ...
                                                        ind_eactive, ...
                                                        tol_g, tol_con);
            end
        end
    else
        break
    end

end


end