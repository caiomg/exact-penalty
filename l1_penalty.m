function [x, history_solution] = l1_penalty(f, phi, x0, mu, epsilon, delta, Lambda)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

global x_prev
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

radius = 1;
radius_max = 100;

history_solution.x = x;
history_solution.rho = NaN;
history_counter = 1;
iter = 0;
finish = false;
[fx, gfx] = f(x);
Hfx = eye(length(x));
[~, ~, Hfx] = f(x);
while ~finish
    iter = iter + 1;
    if iter == 2
        1;
    end
    [ind_eactive, ind_eviolated] = ...
        identify_new_constraints(current_constraints, epsilon, ...
                                 []);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive);
    pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                         ind_eviolated);

    p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);


    if (norm(N'*pseudo_gradient) > max(Lambda, tol_g))
    	x_prev = x;
        B = zeros(size(Hfx));
        for n = ind_eviolated'
            B = B + mu*(current_constraints(n).H);
        end
        B = B + Hfx;
        % Calculate Newton direction
        nBn = N'*B*N;
        if rcond(nBn) > sqrt(eps) && (pseudo_gradient'*(N*N'))*B*((N*N')*pseudo_gradient) > sqrt(eps) 
            u = -(nBn)\(N'*pseudo_gradient);
        else
            u = -(N'*pseudo_gradient);
        end
        h0 = N*u;
        h = correct_direction(h0, Q*R);
        model.B = N'*B*N;
        model.g = N'*pseudo_gradient;
        Ii = zeros(0, 1);
        for constn = 1:n_constraints
            if isempty(find(ind_eactive == constn, 1))
               Ii(end+1, 1) = constn; 
            end
        end
        [s, fs, ind_eactive1] = cauchy_step(model, radius, N, mu, current_constraints, Ii, u, zeros(size(u)), ind_eactive, epsilon);
        p1b = @(x) l1_function(f, phi, mu, x, ind_eactive1);
        pred = fs;
        fmodel.H = Hfx;
        fmodel.g = gfx;
        ared = p(x) - p(x + N*s);
        rho = ared/pred;
        if rho < 0.25
            radius = radius/4;
        else
            if rho > 3/4 && norm(s) > 0.99*radius
                radius = min(2*radius, radius_max);
            end
        end
        if rho > 0.1
            x = x + N*s;
        end
%         x = l1_linear_search(B, pseudo_gradient, h, x, ...
%                               current_constraints, p, mu, epsilon, ...
%                               delta, ind_eactive);
        history_solution(end+1).x = x;
        history_solution(end).rho = rho;
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
            [~, ~, Hfx] = f(x);
        end
        current_constraints = evaluate_constraints(phi, x);
    else
        iter2 = 0;
        while true
            p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);
            pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, ...
                                                 ind_eviolated);

          	x_prev = x;
            [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive);
            % calculate multipliers
            multipliers = -R\(Q'*pseudo_gradient);
    %         multipliers = linsolve(R, Q'*pseudo_gradient, linsolve_opts);

            Bv = zeros(size(Hfx));
            for n = ind_eviolated'
                Bv = Bv + mu*(current_constraints(n).H);
            end
            Bm = zeros(size(Hfx));
            for n = 1:length(ind_qr)
                Bm = Bm + multipliers(n)*(current_constraints(ind_qr(n)).H);
            end
            B = Bv + Bm + Hfx;

            % Are there conditions for dropping one constraint?
            if sum(multipliers < 0 | mu < multipliers)
                [h, sigma, grad_phi_j, Q1, R1, ind_j] = ...
                                    l1_drop_constraint(Q, R, multipliers, mu);
                h = correct_direction(h, Q1*R1);
                ind_null = sum(abs(R1'), 1) < 1e-10;
                N1 = Q1(:, ind_null);
                if ((pseudo_gradient + max(0, sigma)*grad_phi_j)'*h < -delta)
                    n_drop = ind_qr(ind_j);
                    ind_eactive_dropping = ind_eactive(ind_eactive ~= n_drop);
                    p2 = @(x) l1_function(f, phi, mu, x, ind_eactive_dropping);
                    B = B - multipliers(ind_j)*(current_constraints(n_drop).H);
                    if current_constraints(n_drop).c > 0
                        pseudo_gradient = l1_pseudo_gradient(gfx, mu, current_constraints, [ind_eviolated; n_drop]);
                        B = B + mu*(current_constraints(n_drop).H);
                    end
                    model.B = N1'*B*N1;
                    model.g = N1'*pseudo_gradient;
                    Ii = zeros(0, 1);
                    for constn = 1:n_constraints
                        if isempty(find(ind_eactive_dropping == constn, 1))
                           Ii(end+1, 1) = constn; 
                        end
                    end
                    [s, fs, ind_eactive_dropping_b] = cauchy_step(model, radius, N1, mu, current_constraints, Ii, [], zeros(size(model.g)), ind_eactive_dropping, epsilon);
                    p2b = @(x) l1_function(f, phi, mu, x, ind_eactive_dropping_b);
                    pred = fs;
                    model.H = Hfx;
                    model.g = gfx;
                    ared = p(x) - p(x + N1*s);
                    rho = ared/pred;
                    if rho < 0.25
                        radius = radius/4;
                    else
                        if rho > 3/4 && norm(s) > 0.99*radius
                            radius = min(2*radius, radius_max);
                        end
                    end
                    if rho > 0.1
                        x = x + N1*s;
                    end
                    
%                     x = l1_linear_search(B, pseudo_gradient, h, x, ...
%                                   current_constraints, p2, mu, epsilon, ...
%                                   delta, ind_eactive);
                    history_solution(end+1).x = x;
                    history_solution(end).rho = rho;
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
                        [~, ~, Hfx] = f(x);
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
                                        [], tol_g, ...
                                        tol_con);
                end
                break
            elseif (norm(N'*pseudo_gradient) < tol_g && ...
                    isempty(find(multipliers < 0 | multipliers > mu, 1)))
                
                n_qr = size(ind_qr, 1);
                phih = zeros(n_qr, 1);
                for n = 1:n_qr
                   phih(n) = current_constraints(ind_qr(n)).c;
                end
                if norm(phih) < tol_con
                    finish = true;
                    break
                else
%                 % Test complementarity condition
%                 if norm(multipliers.*vertcat(current_constraints(ind_qr).c), 1) > 100*tol_con
                    % Case to lower epsilon and drop constraints
                    [epsilon, Lambda, N, Q, R, ind_eactive, ind_eviolated, ...
                     ind_qr] = l1_criticality_step(epsilon, Lambda, ....
                                                   current_constraints, ...
                                                   gfx, mu, Q, R, ...
                                                   [], tol_g, ...
                                                   tol_con);
                end
%                 else
%                     finish = true;
%                     break
%                 end
            else
                %%%%%%%%%%%%
                % Calculate Newton direction
                if ~isempty(N)
                    if rcond(N'*B*N) > sqrt(eps) && (pseudo_gradient'*(N*N'))*B*((N*N')*pseudo_gradient) > sqrt(eps) 
                        u = -(N'*B*N)\(N'*pseudo_gradient);
                    else
                        u = -(N'*pseudo_gradient);
                    end
                    h = N*u;
                    h = correct_direction(h, Q*R);
                else
                    h = zeros(size(x));
                    u = zeros(0, 1);
                end
                %%%%%%%%%%%%%%%
                model.B = N'*B*N;
                model.g = N'*pseudo_gradient;
                Ii = zeros(0, 1);
                for constn = 1:n_constraints
                    if isempty(find(ind_eactive == constn, 1))
                       Ii(end+1, 1) = constn; 
                    end
                end
                [s, fs, ind_eactive_b] = cauchy_step(model, radius, N, mu, current_constraints, Ii, u, zeros(size(u)), ind_eactive, epsilon);
                p1b = @(x) l1_function(f, phi, mu, x, ind_eactive_b);
                pred1 = fs;
                model.H = Hfx;
                model.g = gfx;
                %%%%%%%%%%%%
                % Recalculate constraints
                [n_qr, ~] = size(ind_qr);
                phih = zeros(n_qr, 1);
                for n = 1:n_qr
                   phih(n) = phi{ind_qr(n)}(x + h);
                end
                %%%%%%%%%%%%%%%%%%                
                v = tr_vertical_step(Q, R, phih, N*s, radius);
                normphi = norm([current_constraints(ind_eactive).c], 1);
                ppgrad = N'*pseudo_gradient;
                % Remaining: compare h+v to h, beforehand
                pred = predict_descent(model, current_constraints, N*s + v, mu, ind_eactive);
                if pred >= delta*(norm(ppgrad)^2 + normphi)
                    ared = p(x) - p(x + N*s + v);
                    rho = ared/pred;
                    if rho < 0.25
                        radius = radius/4;
                    else
                        if rho > 3/4 && norm(s) > 0.99*radius
                            radius = min(2*radius, radius_max);
                        end
                    end
                    if rho > 0.1
                        x = x + N*s + v;
                        current_constraints = evaluate_constraints(phi, x);
                    end
                    history_solution(end+1).x = x;
                    history_solution(end).rho = rho;
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
                        [~, ~, Hfx] = f(x);
                    end
                else
                     [epsilon, Lambda, N, Q, R, ind_eactive, ind_eviolated, ...
                     ind_qr] = l1_criticality_step(epsilon, Lambda, ....
                                                   current_constraints, ...
                                                   gfx, mu, Q, R, ...
                                                   [], tol_g, ...
                                                   tol_con);
                    break
                end

            end
            iter2 = iter2 + 1;
        end
    end
end


end