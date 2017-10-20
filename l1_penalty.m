function x = l1_penalty(f, phi, x0, mu, epsilon, delta, Lambda )
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

linsolve_opts.UT = true;

x = x0;
n_variables = size(x, 1);
n_constraints = size(phi, 1);
tol_g = 1e-6;
tol_con = 1e-6;

p = @(x) l1_function(f, phi, mu, x);
% QR decomposition of constraints gradients matrix A
Q = zeros(n_variables, 0);
R = zeros(0, 0);
ind_eactive = zeros(0, 1);


current_constraints = evaluate_constraints(phi, x);

iter = 0;
finish = false;
[fx, gfx] = f(x);
Hfx = eye(length(x));
while ~finish
    iter = iter + 1;
    [ind_eactive, ind_eviolated] = ...
        identify_new_constraints(current_constraints, epsilon, ...
                                 []);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_eviolated);

    p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);


    if (norm(N'*pseudo_gradient) > Lambda)
    	x_prev = x;
        B = zeros(size(Hfx));
        for n = ind_eviolated'
            B = B - (current_constraints(n).H)/mu;
        end
        B = B + Hfx;
        % Calculate Newton direction
        nBn = N'*B*N;
        if ~issymmetric(nBn)
            nBn = (nBn + nBn')/2;
        end
        if rcond(nBn) > sqrt(eps) && eigs(nBn, 1, 'SA') > sqrt(eps)
            u = -(nBn)\(N'*pseudo_gradient);
        else
            u = -(N'*pseudo_gradient);
        end
        h0 = N*u;
        h = correct_direction(h0, Q*R);
        x = l1_linear_search(B, pseudo_gradient, h, x, ...
                              current_constraints, p, mu, epsilon, ...
                              delta, ind_eactive);
        if norm(x - x_prev, 'inf') ~= 0
            gfx_prev = gfx;
            [fx, gfx] = f(x);
            s = x - x_prev;
            y = gfx - gfx_prev;
            if s'*y >= 0.2*(s'*Hfx*s)
                theta = 1;
            else
                theta = 0.8*(s'*Hfx*s)/(s'*Hfx*s - s'*y);
            end
            Hs = Hfx*s;
            r = theta*y + (1-theta)*(Hs);
            t1 = - ((Hs*Hs')/(s'*Hfx*s));
            t2 = (r*r')/(s'*r);
            Hfx = Hfx + t1 + t2;
        end
        current_constraints = evaluate_constraints(phi, x);
    else
        iter2 = 0;
        while true
          	x_prev = x;
            [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive);
            % calculate multipliers
            multipliers = R\(Q'*pseudo_gradient);
    %         multipliers = linsolve(R, Q'*pseudo_gradient, linsolve_opts);

            B = zeros(size(Hfx));
            for n = ind_eviolated'
                B = B - (current_constraints(n).H)/mu;
            end
            for n = 1:length(ind_qr)
                B = B - multipliers(n)*(current_constraints(ind_qr(n)).H);
            end
            B = B + Hfx;

            % Are there conditions for dropping one constraint?
            if sum(multipliers < 0 | 1/mu < multipliers)
                [h, sigma, grad_phi_j, Q1, R1, ind_j] = ...
                                    l1_drop_constraint(Q, R, multipliers, mu);
                h = correct_direction(h, Q1*R1);
                if ((pseudo_gradient + min(0, sigma)*grad_phi_j)'*h < -delta)
                    n_drop = ind_qr(ind_j);
                    p2 = @(x) l1_function(f, phi, mu, x, ind_eactive(ind_eactive ~= n_drop));
                    B = B + multipliers(ind_j)*(current_constraints(n_drop).H);
                    if current_constraints(n_drop).c < 0
                        B = B - (current_constraints(n_drop).H)/mu;
                    end
                    x = l1_linear_search(B, pseudo_gradient, h, x, ...
                                  current_constraints, p2, mu, epsilon, ...
                                  delta, ind_eactive);
                    if norm(x - x_prev, 'inf') ~= 0
                        gfx_prev = gfx;
                        [fx, gfx] = f(x);
                        s = x - x_prev;
                        y = gfx - gfx_prev;
                        if s'*y >= 0.2*(s'*Hfx*s)
                            theta = 1;
                        else
                            theta = 0.8*(s'*Hfx*s)/(s'*Hfx*s - s'*y);
                        end
                        r = theta*y + (1-theta)*(Hfx*s);
                        Hfx = Hfx - ((Hfx*s)*(s'*Hfx)/(s'*Hfx*s)) + (r*r')/(s'*r);
                    end


                    current_constraints = evaluate_constraints(phi, x);
                    Q = Q1;
                    R = R1;
                    ind_qr = ind_qr(ind_qr ~= n_drop);
                    ind_eactive = ind_eactive(ind_eactive ~= n_drop);
                else
                    [epsilon, Lambda, N, Q, R] = ...
                        l1_criticality_step(epsilon, Lambda, ....
                                        current_constraints, ...
                                        gfx, mu, Q, R, ...
                                        ind_eactive, tol_g, ...
                                        tol_con);
                end
                break
            elseif (norm(N'*pseudo_gradient) < tol_g && ...
                    norm(min(0, vertcat(current_constraints.c)), 1) < tol_con ...
                    && isempty(find(multipliers < 0 | multipliers > 1/mu, 1)))
                % Test complementarity condition
                if norm(multipliers.*vertcat(current_constraints(ind_qr).c), 1) > tol_con
                    % Case to lower epsilon and drop constraints
                    [epsilon, Lambda, N, Q, R, ind_eactive, ind_eviolated, ...
                     ind_qr] = l1_criticality_step(epsilon, Lambda, ....
                                                   current_constraints, ...
                                                   gfx, mu, Q, R, ...
                                                   ind_eactive, tol_g, ...
                                                   tol_con);
                else
                    finish = true;
                    break
                end
            else
                % Calculate Newton direction
                if rcond(N'*B*N) > sqrt(eps) && eigs(N'*B*N, 1, 'SA') > sqrt(eps)
                    u = -(N'*B*N)\(N'*pseudo_gradient);
                else
                    u = -(N'*pseudo_gradient);
                end
                h = N*u;
                h = correct_direction(h, Q*R);
                % Recalculate constraints
                [n_qr, ~] = size(ind_qr);
                phih = zeros(n_qr, 1);
                for n = 1:n_qr
                   phih(n) = phi{ind_qr(n)}(x + h);
                end
                % Vertical step
                v = -Q*(R'\phih);
                normphi = norm([current_constraints(ind_eactive).c], 1);
                ppgrad = N'*pseudo_gradient;
                if (p(x + h + v) <= p(x) - delta*(norm(ppgrad)^2 + normphi))
                    x = x + h + v;
                    current_constraints = evaluate_constraints(phi, x);
                    if norm(x - x_prev, 'inf') ~= 0
                        gfx_prev = gfx;
                        [fx, gfx] = f(x);
                        s = x - x_prev;
                        y = gfx - gfx_prev;
                        if s'*y >= 0.2*(s'*Hfx*s)
                            theta = 1;
                        else
                            theta = 0.8*(s'*Hfx*s)/(s'*Hfx*s - s'*y);
                        end
                        r = theta*y + (1-theta)*(Hfx*s);
                        Hfx = Hfx - ((Hfx*s)*(s'*Hfx)/(s'*Hfx*s)) + (r*r')/(s'*r);
                    end

                else
                    [epsilon, Lambda, N, Q, R, ind_eactive, ind_eviolated, ...
                     ind_qr] = l1_criticality_step(epsilon, Lambda, ....
                                                   current_constraints, ...
                                                   gfx, mu, Q, R, ...
                                                   ind_eactive, tol_g, ...
                                                   tol_con);
                    break
                end
            end
        end
    end
end


end