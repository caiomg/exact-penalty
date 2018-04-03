function [x, history_solution] = l1_penalty_article(f, phi, initial_points, mu, epsilon, delta, Lambda)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

options = struct('tol_radius', 1e-5, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_0', 0, 'eta_1', 0.1, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 1, 'radius_max', 1e3, ...
                        'criticality_mu', 100, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'full quadratic', ...
                        'pivot_threshold', 1/6);

gamma_0 = 0.0625;
gamma_1 = 0.5;
gamma_2 = 2;
ut_option.UT = true;

all_f = {f, phi{:}};
n_functions = size(all_f, 2);
initial_fvalues = [];
[dimension, n_initial_points] = size(initial_points);
% Calculating function values for other points of the set
for nf = 1:n_functions
    for k = 1:n_initial_points
        initial_fvalues(nf, k) = all_f{nf}(initial_points(:, k));
    end
end

% Calculating basis of polynomials
switch options.basis
    case 'linear'
        basis = natural_basis(dimension, dimension+1);
    case 'full quadratic'
        basis = natural_basis(dimension);
    case 'diagonal hessian'
        basis = diagonal_basis(dimension);
end

% Initializing model structure
trmodel.points = initial_points;
trmodel.fvalues = initial_fvalues;
trmodel.radius = options.initial_radius;
trmodel.basis = basis;

fphi = {f, phi{:}}';

% Completing set of interpolation points and calculating polynomial
% model
% trmodel = complete_interpolation_set(trmodel, fphi, options);

fval_current = trmodel.fvalues(1, 1);


global x_prev
linsolve_opts.UT = true;

x0 = initial_points(:, 1);
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

radius_max = 100;
history_solution.x = x;
rho = nan;
history_solution.rho = rho;
history_solution.radius = trmodel.radius;

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
                                                  Q, R, ind_eactive, false);
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                         ind_eviolated);

    p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);

    Bv = zeros(dimension);
    for n = ind_eviolated'
        Bv = Bv + mu*(current_constraints(n).H);
    end
    B = Bv + fmodel.H;
    if (norm(N'*pseudo_gradient) > max(Lambda, tol_g))
        x_prev = x;


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
            [d, fd] = solve_tr_problem(model.B, model.g, rf*trmodel.radius);
            [s, fs, ind_eactive1] = cauchy_step(model, trmodel.radius, N, ...
                                                mu, current_constraints, ...
                                                Ii, d, zeros(size(u)), ...
                                                ind_eactive, ...
                                                epsilon);
            
            Ns = N*s;
            pred_r = predict_descent(fmodel, current_constraints, Ns, mu, ind_qr);
            
            if pred_r > 0% && norm(s) > 0.1*norm(d)
                break
            else
                gamma_dec = 0.75;%max(0.5, gamma_1*norm(Ns)/(rf*trmodel.radius));
                rf = gamma_dec*rf;
                if rf < 1e-100
                    1;
                end
            end
        end

%         v = zeros(size(v));
        v = tr_vertical_step_new(fmodel, current_constraints, mu, Ns, ind_eactive, ind_eviolated, trmodel.radius);
        normphi = norm([current_constraints(ind_eactive).c], 1);
        ppgrad = N'*pseudo_gradient;
        pred = predict_descent(fmodel, current_constraints, Ns + v, mu, []);
        if pred < 0
            pred = pred_r;
        end
%%%%%%%%%%%%%
        step = Ns + v;
%         [step, pred] = line_search_full_domain(fmodel, current_constraints, mu, Ns + v, trmodel.radius);
        trial_point = x + step;
        p_trial = p(trial_point);
        ared = px - p_trial;
        if ared < 0
            dpred = pred - 10*eps*max(1, abs(px));
            dared = ared - 10*eps*max(1, abs(px));
        else
            dpred = pred;
            dared = ared;
        end
        if abs(dared) < 10*eps && abs(dpred) < 10*eps
            rho = 1;
        else
            rho = dared/dpred;
        end
        if rho < 0.25
            gamma_dec = max(gamma_0, gamma_1*norm(step)/trmodel.radius);
            trmodel.radius = gamma_dec*trmodel.radius;
        else
            if rho > 3/4
                radius_inc = max(1, gamma_2*(norm(step)/trmodel.radius));
                trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
            end
        end
        if rho > 0.1
            x = trial_point;
            step_accepted = true;
%             trmodel = move_trust_region(trmodel, x, new_center_fvals, ...
%                                    fphi, options)
            px = p_trial;
            [~, fmodel.g, fmodel.H] = f(x);
            current_constraints = evaluate_constraints(phi, x);
        else
            step_accepted = false;
%             trmodel = try_to_add_interpolation_point(trmodel, x, ...
%                                                 new_point_fvals, ...
%                                                 fphi, options)
        end
%         if abs(rho - 1) < 1e-4 && norm(step) < 0.1*norm(N*s)
%             Lambda = Lambda*1.5;
%         end
        history_solution(end+1).x = x;
        history_solution(end).rho = rho;
        history_solution(end).radius = trmodel.radius;
    else
        iter2 = 0;
        while true
            p1 = @(x) l1_function(f, phi, mu, x, ind_eactive);
            pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                                 ind_eviolated);

          	x_prev = x;
            [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive, true);
            % calculate multipliers
            rows_qr = size(R, 1) - size(N, 2);
            multipliers_a = -linsolve(R(1:rows_qr, :), (Q(:, 1:rows_qr)'*pseudo_gradient), ut_option);
            remain = -((Q*R)*multipliers_a + pseudo_gradient);
            correction = linsolve(R(1:rows_qr, :), Q(:, 1:rows_qr)'*remain, ut_option);
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
                [d, fd] = solve_tr_problem(model.B, model.g, trmodel.radius);
                [s, fs, ind_eactive_dropping_b] = cauchy_step(model, trmodel.radius, N1, mu, current_constraints, Ii, d, zeros(size(model.g)), ind_eactive_dropping, epsilon);
                Ns = N1*s;
                pred1 = predict_descent_with_multipliers(fmodel, ...
                                                         current_constraints, Ns, mu, ind_qr_dropping, multipliers_dropping);
%                 [step, pred] = line_search_full_domain(fmodel, ...
%                                                        current_constraints, ...
%                                                        mu, Ns, trmodel.radius);
                pred = pred1;
                step = Ns;
                if pred > delta
                    dropping_succeeded = true;
%                     step = Ns;
                    trial_point = x + step;
                    p_trial = p(trial_point);
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
                        gamma_dec = max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        if rho > 3/4
                            radius_inc = max(1, gamma_2*(norm(step)/trmodel.radius));
                            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                        end
                    end
                    if rho > 0.1
                        x = trial_point;
                        step_accepted = true;
                        px = p_trial;
                        [~, fmodel.g, fmodel.H] = f(x);
                        current_constraints = evaluate_constraints(phi, x);
                    else
                        step_accepted = false;
                    end
                    Q = Q1;
                    R = R1;
                    ind_qr = ind_qr(ind_qr ~= n_drop);
                    ind_eactive = ind_eactive(ind_eactive ~= n_drop);
                else % min decrease not satisfied
                    dropping_succeeded = false;
                    step_accepted = false;
                    gamma_dec = max(gamma_0, gamma_1*norm(s)/trmodel.radius);
                    trmodel.radius = gamma_dec*trmodel.radius;
                end
                history_solution(end+1).x = x;
                history_solution(end).rho = rho;
                history_solution(end).radius = trmodel.radius;


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
                model.B = N'*B*N;
                model.g = N'*pseudo_gradient;
                Ii = zeros(0, 1);
                for constn = 1:n_constraints
                    if isempty(find(ind_eactive == constn, 1))
                       Ii(end+1, 1) = constn; 
                    end
                end
                step_calculation_ok = true;
                d = solve_tr_problem(model.B, model.g, trmodel.radius);
                [s, fs, ind_eactive_b] = cauchy_step(model, ...
                                                     trmodel.radius, ...
                                                     N, mu, ...
                                                     current_constraints, ...
                                                     Ii, d, ...
                                                     zeros(size(model.g)), ...
                                                     ind_eactive, epsilon);
                pred_h = predict_descent(fmodel, current_constraints, N*s, mu, []);
                v = tr_vertical_step_new(fmodel, current_constraints, ...
                                         mu, N*s, ind_eactive, ...
                                         ind_eviolated, trmodel.radius);
                normphi = norm([current_constraints(ind_eactive).c], 1);
                ppgrad = N'*pseudo_gradient;
                pred = predict_descent(fmodel, current_constraints, ...
                                       N*s + v, mu, []);
%                 if pred < delta*(norm(ppgrad)^2 + normphi)
% %                     [Ns, pred] = line_search_full_domain(fmodel, current_constraints, mu, N*s, trmodel.radius);
%                     v = tr_vertical_step_new(fmodel, current_constraints, ...
%                                              mu, Ns, ind_eactive, ...
%                                              ind_eviolated, trmodel.radius);
%                     pred = predict_descent(fmodel, current_constraints, ...
%                                            Ns + v, mu, []);
%                 end
                if pred < delta*(norm(ppgrad)^2 + normphi)
                    % Better not to try now
                    step_accepted = false;
                    gamma_dec = max(gamma_0, gamma_1*norm(s)/trmodel.radius);
                    trmodel.radius = gamma_dec*trmodel.radius;
                else
                    % Compute ared and all...
                    p2 = @(x) l1_function_2nd_order(f, phi, mu, x, [], multipliers, ind_qr);
                    step = N*s + v;
                    trial_point = x + step;
                    p_trial = p(trial_point);
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
                        gamma_dec = max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        if rho > 3/4
                            radius_inc = max(1, gamma_2*(norm(step)/trmodel.radius));
                            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                        end
                    end
                    if rho > 0.1
                        x = trial_point;
                        current_constraints = evaluate_constraints(phi, x);
                        step_accepted = true;
                        px = p_trial;
                        [~, fmodel.g, fmodel.H] = f(x);
                    else
                        step_accepted = false;
                    end
                end
                history_solution(end+1).x = x;
                history_solution(end).rho = rho;
                history_solution(end).radius = trmodel.radius;

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