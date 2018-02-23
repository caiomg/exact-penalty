function [x, history_solution] = l1_penalty_article(f, phi, x0, mu, epsilon, delta, Lambda)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

global x_prev
linsolve_opts.UT = true;

x = x0;
dimension = size(x, 1);
n_constraints = size(phi, 1);
tol_g = 1e-6;
tol_con = 1e-6;
gamma1 = 0.01;

p = @(x) l1_function(f, phi, mu, x);
% QR decomposition of constraints gradients matrix A
Q = zeros(dimension, 0);
R = zeros(0, 0);
ind_eactive = zeros(0, 1);


current_constraints = evaluate_constraints(phi, x);

radius = 1;
radius_max = 100;

history_solution.x = x;
history_solution.rho = NaN;
iter = 0;
finish = false;
[~, fmodel.g] = f(x);
fmodel.H = eye(length(x));
[~, ~, fmodel.H] = f(x);

fmodel.g = fmodel.g;
fmodel.H = fmodel.H;
px = p(x);
while ~finish
    iter = iter + 1;

    [ind_eactive, ind_eviolated] = ...
        identify_new_constraints(current_constraints, epsilon, ...
                                 []);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive);
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                         ind_eviolated);

    p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);


    if (norm(N'*pseudo_gradient) > max(Lambda, tol_g))
        x_prev = x;
        B = zeros(dimension);
        for n = ind_eviolated'
            B = B + mu*(current_constraints(n).H);
        end
        B = B + fmodel.H;
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

        step_calculation_ok = true;
        rf = 1;
        while step_calculation_ok
            [s, fs, ind_eactive1] = cauchy_step(model, rf*radius, N, ...
                                                mu, current_constraints, ...
                                                Ii, [], zeros(size(u)), ...
                                                ind_eactive, ...
                                                epsilon);
            Ns = N*s;
            Ns = correct_direction(Ns, Q*R);
            pv = @(s) -predict_descent(fmodel, current_constraints, s, mu);
            Ns = simple_backtracking_line_search(pv, zeros(size(x)), Ns, 0, 0.9, 50);

            %%%%%%%%%%%%
            % Predict constraint values
            [n_qr, ~] = size(ind_qr);
            phih = zeros(n_qr, 1);
            for n = 1:n_qr
               phih(n) = current_constraints(ind_qr(n)).c + current_constraints(ind_qr(n)).g'*Ns + 0.5*((s'*N')*current_constraints(ind_qr(n)).H*(Ns));
            end
            %%%%%%%%%%%%%%%%%%
            pv = @(s) -predict_descent(fmodel, current_constraints, s, mu);
            v1 = tr_vertical_step(pv, x, Q, R, phih, Ns, radius);
%             v2 = tr_new_vertical_step(pv, current_constraints, Ns, rf*radius, ind_qr, fmodel);
%             v3 = tr_new_vertical_step_alternative(pv, current_constraints, Ns, rf*radius, ind_qr, fmodel, mu);
            v = v1;
            normphi = norm([current_constraints(ind_eactive).c], 1);
            ppgrad = N'*pseudo_gradient;
            pred = predict_descent(fmodel, current_constraints, Ns + v, mu, []);
            if pred > 0
                break
            else
                rf = min(rf/2, norm(Ns)/norm(s));
            end
        end
%%%%%%%%%%%%%
        p_trial = p(x + Ns + v);
        ared = px - p_trial;
        dpred = pred - 10*eps*max(1, abs(px));
        dared = ared - 10*eps*max(1, abs(px));
        if abs(dared) < 10*eps && abs(dpred) < 10*eps
            rho = 1;
        else
            rho = dared/dpred;
        end
        if rho < 0.25
            radius = radius/4;
        else
            if rho > 3/4 && norm(s) > 0.85*radius
                radius = min(2*radius, radius_max);
            end
        end
        if rho > 0.1
            x = x + Ns + v;
            step_accepted = true;
            px = p_trial;
            [~, fmodel.g, fmodel.H] = f(x);
            current_constraints = evaluate_constraints(phi, x);
        else
            step_accepted = false;
        end
        history_solution(end+1).x = x;
        history_solution(end).rho = rho;
    else
        iter2 = 0;
        while true
            p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);
            pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                                 ind_eviolated);

          	x_prev = x;
            [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive);
            % calculate multipliers
            multipliers_a = -R\(Q'*pseudo_gradient);
            remain = -(Q*R*multipliers_a + pseudo_gradient);
            correction = R\(Q'*remain);
            multipliers = multipliers_a + correction;
            tol_multipliers = 10*norm(correction);
    %         multipliers = linsolve(R, Q'*pseudo_gradient, linsolve_opts);

            Bv = zeros(dimension);
            for n = ind_eviolated'
                Bv = Bv + mu*(current_constraints(n).H);
            end
            Bm = zeros(dimension);
            for n = 1:length(ind_qr)
                Bm = Bm + multipliers(n)*(current_constraints(ind_qr(n)).H);
            end
            B = Bv + Bm + fmodel.H;

            dropping_constraint = false;
            % Are there conditions for dropping one constraint?
            if sum(multipliers < -tol_multipliers | mu < multipliers - tol_multipliers)
                dropping_constraint = true;
                [h, sigma, grad_phi_j, Q1, R1, ind_j] = ...
                                    l1_drop_constraint(Q, R, multipliers, mu);
                h = correct_direction(h, Q1*R1);
                ind_null = sum(abs(R1'), 1) < 1e-10;
                N1 = Q1(:, ind_null);
                n_drop = ind_qr(ind_j);
                ind_eactive_dropping = ind_eactive(ind_eactive ~= n_drop);
                ind_qr_dropping = ind_qr(ind_qr ~= n_drop);
                multipliers_dropping = multipliers([(1:ind_j-1)';(ind_j+1:end)']);
                p2 = @(x) l1_function(f, phi, mu, x, ind_eactive_dropping);
                B = B - multipliers(ind_j)*(current_constraints(n_drop).H);
                if current_constraints(n_drop).c > 0
                    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, [ind_eviolated; n_drop]);
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
                contador = 0;
                step_calculation_ok = true;
                rf = 1;
                while step_calculation_ok
                    [s, fs, ind_eactive_dropping_b] = cauchy_step(model, rf*radius, N1, mu, current_constraints, Ii, [], zeros(size(model.g)), ind_eactive_dropping, epsilon);
                    Ns = N1*s;
                    Ns = correct_direction(Ns, Q1*R1);
                    pv = @(s) -predict_descent_with_multipliers(fmodel, current_constraints, s, mu, ind_qr_dropping, multipliers_dropping);
                    Ns = simple_backtracking_line_search(pv, zeros(size(x)), Ns, 0, 0.9, 50);
                    pred1 = predict_descent_with_multipliers(fmodel, current_constraints, Ns, mu, ind_qr_dropping, multipliers_dropping);
                    pred = predict_descent(fmodel, current_constraints, Ns, mu, []);
                    if pred > 0
                        break
                    else
                        rf = min(rf/2, norm(Ns)/norm(s));
                        contador = contador + 1;
                    end
                    if norm(Ns)/norm(s) < 0.1
                        step_calculation_ok = false;
                    end
                end
                if step_calculation_ok
                    p2 = @(x) l1_function_2nd_order(f, phi, mu, x, [], multipliers_dropping, ind_qr_dropping);
                    p_trial = p(x + Ns);
                    ared = px - p_trial;
                    ared1 = p2(x) - p2(x + Ns);
                    dpred = pred - 10*eps*max(1, abs(px));
                    dared = ared - 10*eps*max(1, abs(px));
                    if abs(dared) < 10*eps && abs(dpred) < 10*eps
                        rho = 1;
                    else
                        rho = dared/dpred;
                    end
                    if rho < 0.25
                        radius = radius/4;
                    else
                        if rho > 3/4 && norm(s) > 0.85*radius
                            radius = min(2*radius, radius_max);
                        end
                    end
                    if rho > 0.1
                        x = x + Ns;
                        step_accepted = true;
                        px = p_trial;
                        [~, fmodel.g, fmodel.H] = f(x);
                    else
                        step_accepted = false;
                    end
                    % Testing 'line-search' condition
                    if (N*(N'*current_constraints(n_drop).g))'*pseudo_gradient > delta
                        dropping_succeeded = true;
                    else
                        dropping_succeeded = false;
                    end
                    if step_accepted
                        current_constraints = evaluate_constraints(phi, x);
                    end
                else
                    dropping_succeeded = false;
                    % The l1 criticality step will be used
                end

                history_solution(end+1).x = x;
                history_solution(end).rho = rho;

                Q = Q1;
                R = R1;
                ind_qr = ind_qr(ind_qr ~= n_drop);
                ind_eactive = ind_eactive(ind_eactive ~= n_drop);
            elseif (norm(N'*pseudo_gradient) < tol_g)% && ...
                    %~isempty(find(multipliers > -10*eps & multipliers < mu + 10*(eps(mu)), 1)))
                n_qr = size(ind_qr, 1);
                phih = zeros(n_qr, 1);
                for n = 1:n_qr
                   phih(n) = current_constraints(ind_qr(n)).c;
                end
                if norm(phih) < tol_con
                    finish = true;
                    break
                else
                    step_accepted = false;
                end
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
                [s, fs, ind_eactive_b] = cauchy_step(model, radius, N, mu, current_constraints, Ii, [], zeros(size(u)), ind_eactive, epsilon);
                p1b = @(x) l1_function(f, phi, mu, x, ind_eactive_b);
                pred_h = fs;
                pred_h = predict_descent(fmodel, current_constraints, N*s, mu, []);
                %%%%%%%%%%%%
                % Predict constraint values
                [n_qr, ~] = size(ind_qr);
                phih = zeros(n_qr, 1);
                for n = 1:n_qr
                   phih(n) = current_constraints(ind_qr(n)).c + current_constraints(ind_qr(n)).g'*N*s + 0.5*((s'*N')*current_constraints(ind_qr(n)).H*(N*s));
                end
                %%%%%%%%%%%%%%%%%%                
                
                v = tr_vertical_step(p, x, Q, R, phih, N*s, radius);
                pv = @(s) -predict_descent(fmodel, current_constraints, s, mu);
                v = tr_new_vertical_step(pv, current_constraints, N*s, radius, ind_qr);
                normphi = norm([current_constraints(ind_eactive).c], 1);
                ppgrad = N'*pseudo_gradient;
                pred_hv = predict_descent(fmodel, current_constraints, N*s + v, mu, []);
                if pred_hv < delta*(norm(ppgrad)^2 + normphi) && pred_hv < pred_h
                    pred = pred_h;
                    v = zeros(size(v));
                else
                    pred = pred_hv;
                end
                p2 = @(x) l1_function_2nd_order(f, phi, mu, x, [], multipliers, ind_qr);
                p_trial = p(x + N*s + v);
                ared = px - p_trial;
                ared1 = p2(x) - p2(x + N*s + v);
                dpred = pred - 10*eps*max(1, abs(px));
                dared = ared - 10*eps*max(1, abs(px));
                if abs(dared) < 10*eps && abs(dpred) < 10*eps
                    rho = 1;
                else
                    rho = dared/dpred;
                end
                if rho < 0.25
                    radius = radius/4;
                else
                    if rho > 3/4 && norm(s) > 0.85*radius
                        radius = min(2*radius, radius_max);
                    end
                end
                if rho > 0.1
                    x = x + N*s + v;
                    current_constraints = evaluate_constraints(phi, x);
                    step_accepted = true;
                    px = p_trial;
                    [~, fmodel.g, fmodel.H] = f(x);
                else
                    step_accepted = false;
                end
                history_solution(end+1).x = x;
                history_solution(end).rho = rho;

            end
            iter2 = iter2 + 1;
            if ~step_accepted || (dropping_constraint && ~dropping_succeeded)
                [epsilon, Lambda, N, Q, R, ind_eactive, ind_eviolated, ...
                 ind_qr] = l1_criticality_step(epsilon, Lambda, ....
                                               current_constraints, ...
                                               fmodel.g, mu, Q, R, ...
                                               [], tol_g, ...
                                               tol_con); 
               break
            end
            if dropping_constraint
                break
            end
        end
    end

end


end