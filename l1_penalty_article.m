function [x, history_solution] = l1_penalty_article(f, phi, initial_points, mu, epsilon, delta, Lambda)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

options = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_1', 0, 'eta_2', 0.1, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 0.5, 'radius_max', 1e4, ...
                        'criticality_mu', 50, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'full quadratic', ...
                        'pivot_threshold', 0.2);

gamma_0 = 0.0625;
gamma_1 = 0.5;
gamma_2 = 2;
ut_option.UT = true;
eta_1 = options.eta_1;
eta_2 = options.eta_2;

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

options.basis = 'diagonal hessian';
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
trmodel = complete_interpolation_set(trmodel, fphi, options);

fval_current = trmodel.fvalues(1, 1);


global x_prev
linsolve_opts.UT = true;

x0 = initial_points(:, 1);
x = x0;
dimension = size(x, 1);
n_constraints = size(phi, 1);
tol_g = 1e-5;
tol_con = 1e-5;
tol_radius = options.tol_radius;
gamma1 = 0.01;

p = @(x) l1_function(f, phi, mu, x);
% QR decomposition of constraints gradients matrix A
Q = zeros(dimension, 0);
R = zeros(0, 0);
ind_eactive = zeros(0, 1);



radius_max = options.radius_max;

iter = 0;
finish = false;
[fx, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
current_constraints = extract_constraints_from_tr_model(trmodel);
px = fx + mu*sum(max(0, [current_constraints.c]));
history_solution.x = x;
rho = nan;
history_solution.rho = rho;
history_solution.radius = trmodel.radius;
history_solution.px = px;
history_solution.fx = fx;
while ~finish
    [~, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
    current_constraints = extract_constraints_from_tr_model(trmodel);

    [ind_eactive, ind_eviolated] = ...
        identify_new_constraints(current_constraints, epsilon, ...
                                 []);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive, false);
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                         ind_eviolated);
    ppgrad = N*(N'*pseudo_gradient);
    Bv = zeros(dimension);
    for n = ind_eviolated'
        Bv = Bv + mu*(current_constraints(n).H);
    end
    Ba = zeros(dimension);
    for n = ind_qr'
        if ppgrad'*current_constraints(n).H*ppgrad > 0
            Ba = Ba + mu*current_constraints(n).H;
        end
    end
    B = (Bv + fmodel.H);
    if (norm(N'*pseudo_gradient) > max(Lambda, tol_g))
        x_prev = x;

        model.B = N'*(B + Ba)*N;
        model.g = N'*pseudo_gradient;
%         if ~isempty(ind_qr)
%             A2 = [current_constraints(ind_qr).g];
%             g2 = correct_direction(pseudo_gradient, A2);
%             model.g = N'*g2;
%         end

        Ii = zeros(0, 1);
        for constn = 1:n_constraints
            if isempty(find(ind_eactive == constn, 1))
               Ii(end+1, 1) = constn; 
            end
        end

        geometry_ok = is_lambda_poised(trmodel, options);
        d1 = -model.g*(trmodel.radius/norm(model.g));

        [h1, pred, status] = line_search_full_domain(fmodel, current_constraints, mu, N*d1, trmodel.radius);
        if ~status
            h1 = zeros(size(h1));
        end
        v1 = tr_vertical_step_new(fmodel, current_constraints, mu, h1, ind_eactive, ind_eviolated, trmodel.radius);
        v1 = zeros(size(v1));
%         h1 = zeros(size(h1));
        h1_fmodel = shift_model(fmodel, h1 + v1);
        for m = 1:n_constraints
            h1_constraints(m) = shift_model(current_constraints(m), h1 + v1);
        end
        step = h1 + v1;
        if norm(h1 + v1) < trmodel.radius - tol_radius
            d = solve_tr_problem(model.B, model.g, trmodel.radius);
            s2 = cauchy_step(model, trmodel.radius, N, ...
                                            mu, current_constraints, ...
                                            Ii, d, zeros(size(d)), ...
                                            ind_eactive, ...
                                            epsilon);
            h2 = N*s2;
            v2 = tr_vertical_step_new(fmodel, current_constraints, mu, h2, ind_eactive, ind_eviolated, trmodel.radius);
%             v2 = zeros(size(v2));
            [hv, ~, status] = line_search_full_domain(h1_fmodel, h1_constraints, mu, (h2 + v2) - (h1 + v1), trmodel.radius - norm(h1 + v1));
            if status
               step = (h1 + v1) + hv;
            end
        end
        pred = predict_descent(fmodel, current_constraints, step, mu, []);
%%%%%%%%%%%%%
        if pred < 0 || norm(step) < 0.0625*trmodel.radius
            rho = -inf;
            if geometry_ok
                trmodel.radius = 0.125*trmodel.radius;
            else
                trmodel = improve_model(trmodel, fphi, options);
            end
        else
            trial_point = x + step;
            [p_trial, trial_fvalues] = p(trial_point);
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

            if rho > eta_2 || (rho > eta_1 && geometry_ok)
                x = trial_point;
                step_accepted = true;
                trmodel = move_trust_region(trmodel, x, trial_fvalues, ...
                                       fphi, options);
                px = p_trial;
            else
                step_accepted = false;
                trmodel = try_to_add_interpolation_point(trmodel, trial_point, ...
                                                    trial_fvalues, ...
                                                    fphi, options);
            end
            if rho < eta_2 && geometry_ok
                gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                trmodel.radius = gamma_dec*trmodel.radius;
            else
                if rho > eta_2
                    radius_inc = max(1, gamma_2*(norm(step)/trmodel.radius));
                    trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                end
            end
        end
    else
        while true
          	x_prev = x;
            [~, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
            current_constraints = extract_constraints_from_tr_model(trmodel);
            pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                                 ind_eviolated);

            [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive, true);
            % calculate multipliers
            rows_qr = size(R, 1) - size(N, 2);
            multipliers_a = -linsolve(R(1:rows_qr, :), (Q(:, 1:rows_qr)'*pseudo_gradient), ut_option);
            remain = -((Q*R)*multipliers_a + pseudo_gradient);
            correction = linsolve(R(1:rows_qr, :), Q(:, 1:rows_qr)'*remain, ut_option);
            multipliers = multipliers_a + correction;
            tol_multipliers = 10*max(norm(correction), eps);
    %         multipliers = linsolve(R, Q'*pseudo_gradient, linsolve_opts);

            if (norm(N'*pseudo_gradient) < 100*tol_g)
                [trmodel, epsilon] = tr_criticality_step(trmodel, fphi, epsilon, ...
                                              mu, options);
                [~, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
                current_constraints = ...
                    extract_constraints_from_tr_model(trmodel);
                [ind_eactive, ind_eviolated] = ...
                    identify_new_constraints(current_constraints, epsilon, []);
                [N, Q, R, ind_qr] = update_factorization(current_constraints, [], ...
                                                         [], ...
                                                         ind_eactive, true);
                pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, ...
                                                     current_constraints, ...
                                                     ind_eviolated);
                rows_qr = size(R, 1) - size(N, 2);
                multipliers_a = -linsolve(R(1:rows_qr, :), (Q(:, 1: ...
                                                              rows_qr)'*pseudo_gradient), ut_option);
                remain = -((Q*R)*multipliers_a + pseudo_gradient);
                correction = linsolve(R(1:rows_qr, :), Q(:, 1: ...
                                                         rows_qr)'*remain, ...
                                      ut_option);
                multipliers = multipliers_a + correction;
                tol_multipliers = 10*max(norm(correction), eps);
                
                if ~sum(multipliers < -tol_multipliers |...
                        mu < multipliers - tol_multipliers) 
                    phih = [current_constraints(ind_qr).c]';
                    if norm(N'*pseudo_gradient) < tol_g
                        if norm(phih) < tol_con
                            finish = true;
                            break
                        end
                    end
                end
            end
            
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
                multipliers_dropping = multipliers;
                ind_qr_dropping = ind_qr;
                ind_eactive_dropping = ind_eactive;
                Q1 = Q;
                R1 = R;
                while sum(multipliers_dropping < -tol_multipliers | mu < multipliers_dropping - tol_multipliers)
                    %%%
                    dropping_constraint = true;
                    [h, sigma, grad_phi_j, Q1, R1, ind_j] = ...
                                        l1_drop_constraint(Q1, R1, multipliers_dropping, mu, tol_multipliers);
                    ind_null = sum(abs(R1'), 1) < 1e-10;
                    N1 = Q1(:, ind_null);
                    n_drop = ind_qr_dropping(ind_j);
                    ind_eactive_dropping = ind_eactive_dropping(ind_eactive_dropping ~= n_drop);
                    ind_qr_dropping = ind_qr_dropping(ind_qr_dropping ~= n_drop);
                    %%%%
                    rows_qr = size(R1, 1) - size(N1, 2);
                    multipliers_dropping = -linsolve(R1(1:rows_qr, :), (Q1(:, 1:rows_qr)'*pseudo_gradient), ut_option);
                    remain = -((Q1*R1)*multipliers_dropping + pseudo_gradient);
                    correction = linsolve(R1(1:rows_qr, :), Q1(:, 1:rows_qr)'*remain, ut_option);
                    multipliers_dropping = multipliers_dropping + correction;
                    tol_multipliers = 10*max(norm(correction), eps);
                    %%%%
                    B = B - multipliers(ind_j)*(current_constraints(n_drop).H);
                    if current_constraints(n_drop).c > 0
                        pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, [ind_eviolated; n_drop]);
                        B = B + mu*(current_constraints(n_drop).H);
                    end
                end
                model.B = N1'*B*N1;
                model.g = N1'*pseudo_gradient;
                Ii = zeros(0, 1);
                for constn = 1:n_constraints
                    if isempty(find(ind_eactive_dropping == constn, 1))
                       Ii(end+1, 1) = constn; 
                    end
                end
                geometry_ok = is_lambda_poised(trmodel, options);
                d = solve_tr_problem(model.B, model.g, trmodel.radius);
                % d = truncated_cg_step(model.B, model.g, trmodel.radius);
                % [s, fs, ind_eactive_dropping_b] = cauchy_step(model, trmodel.radius, N1, mu, current_constraints, Ii, d, zeros(size(model.g)), ind_eactive_dropping, epsilon);
                % Ns = N1*s;
                [Ns, pred, status] = line_search_full_domain(fmodel, current_constraints, mu, N1*d, trmodel.radius);
                pred = predict_descent(fmodel, current_constraints, Ns, mu, []);
                step = Ns;
                if pred > delta
                    dropping_succeeded = true;
                    trial_point = x + step;
                    [p_trial, trial_fvalues] = p(trial_point);
                    ared = px - p_trial;
                    dpred = pred - 10*eps*max(1, abs(px));
                    dared = ared - 10*eps*max(1, abs(px));
                    if abs(dared) < 10*eps && abs(dpred) < 10*eps
                        rho = 1;
                    else
                        rho = dared/dpred;
                    end
                    if rho > eta_2 || (rho > eta_1 && geometry_ok)
                        x = trial_point;
                        step_accepted = true;
                        px = p_trial;
                        trmodel = move_trust_region(trmodel, x, ...
                                                    trial_fvalues, ...
                                                    fphi, options);
                    else
                        step_accepted = false;
                        trmodel = try_to_add_interpolation_point(trmodel, ...
                                                                 trial_point, ...
                                                                 trial_fvalues, fphi, options);
                    end
                    if rho < eta_2 && geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        if rho > eta_2
                            radius_inc = max(1, gamma_2*(norm(step)/trmodel.radius));
                            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                        end
                    end
                    Q = Q1;
                    R = R1;
                    ind_qr = ind_qr(ind_qr ~= n_drop);
                    ind_eactive = ind_eactive(ind_eactive ~= n_drop);
                else % min decrease not satisfied
                    dropping_succeeded = false;
                    step_accepted = false;
                    if geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(s)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        trmodel = improve_model(trmodel, fphi, options);
                    end
                end
            else
                geometry_ok = is_lambda_poised(trmodel, options);
                if ~isempty(N)
                    model.B = N'*B*N;
                    model.g = N'*pseudo_gradient;
                    Ii = zeros(0, 1);
                    for constn = 1:n_constraints
                        if isempty(find(ind_eactive == constn, 1))
                           Ii(end+1, 1) = constn; 
                        end
                    end
%%%%%%%%%%%%
%                     d = solve_tr_problem(model.B, model.g, trmodel.radius);
%                     % d = truncated_cg_step(model.B, model.g, trmodel.radius);
%                     [s, fs, ind_eactive_b] = cauchy_step(model, ...
%                                                          trmodel.radius, ...
%                                                          N, mu, ...
%                                                          current_constraints, ...
%                                                          Ii, d, ...
%                                                          zeros(size(model.g)), ...
%                                                          ind_eactive, epsilon);
%                     Ns = N*s;
%                 else
%                     Ns = zeros(dimension, 1);
%                 end
%                 v = tr_vertical_step_new(fmodel, current_constraints, ...
%                                          mu, Ns, ind_eactive, ...
%                                          ind_eviolated, trmodel.radius);
%                 pred = predict_descent(fmodel, current_constraints, ...
%                                        Ns + v, mu, []);
%                 step = Ns + v;
%%%%%%%%%%%%%%%
                    d1 = -model.g*(trmodel.radius/norm(model.g));
                    [h1, pred, status] = line_search_full_domain(fmodel, current_constraints, mu, N*d1, trmodel.radius);
                    if ~status
                        h1 = zeros(size(h1));
                    end
                    v1 = tr_vertical_step_new(fmodel, current_constraints, mu, h1, ind_eactive, ind_eviolated, trmodel.radius);
                    h1_fmodel = shift_model(fmodel, h1 + v1);
                    for m = 1:n_constraints
                        h1_constraints(m) = shift_model(current_constraints(m), h1 + v1);
                    end
                    step = h1 + v1;
                    if norm(h1 + v1) < trmodel.radius - tol_radius
                        d = solve_tr_problem(model.B, model.g, trmodel.radius);
                        s2 = cauchy_step(model, trmodel.radius, N, ...
                                                        mu, current_constraints, ...
                                                        Ii, d, zeros(size(d)), ...
                                                        ind_eactive, ...
                                                        epsilon);
                        h2 = N*s2;
                        v2 = tr_vertical_step_new(fmodel, current_constraints, mu, h2, ind_eactive, ind_eviolated, trmodel.radius);
                        [hv, ~, status] = line_search_full_domain(h1_fmodel, h1_constraints, mu, (h2 + v2) - (h1 + v1), trmodel.radius - norm(h1 + v1));
                        if status
                           step = (h1 + v1) + hv;
                        end
                    end
    %%%%%%%%%%%%%%%%
                else
                    v = tr_vertical_step_new(fmodel, current_constraints, mu, zeros(dimension, 1), ind_eactive, ind_eviolated, trmodel.radius);
                    step = v;
                end
                pred = predict_descent(fmodel, current_constraints, step, mu, []);
                normphi = norm([current_constraints(ind_eactive).c], 1);
                ppgrad = N'*pseudo_gradient;
                if pred < delta*(norm(ppgrad)^2 + normphi)
                    % Better not to try now
                    step_accepted = false;
                    if geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(Ns+v)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        trmodel = improve_model(trmodel, fphi, options);
                    end
                else
                    % Compute ared and all...
                    trial_point = x + step;
                    [p_trial, trial_fvalues] = p(trial_point);
                    ared = px - p_trial;
                    dpred = pred - 10*eps*max(1, abs(px));
                    dared = ared - 10*eps*max(1, abs(px));
                    if abs(dared) < 10*eps && abs(dpred) < 10*eps
                        rho = 1;
                    else
                        rho = dared/dpred;
                    end
                    if rho > eta_2 || (rho > eta_1 && geometry_ok)
                        x = trial_point;
                        step_accepted = true;
                        px = p_trial;
                        trmodel = move_trust_region(trmodel, x, ...
                                                    trial_fvalues, fphi, options);
                    else
                        step_accepted = false;
                        trmodel = try_to_add_interpolation_point(trmodel, trial_point, ...
                                                                 trial_fvalues, ...
                                                                 fphi, options);
                    end
                    if rho < eta_2 && geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        if rho > eta_2
                            radius_inc = max(1, gamma_2*(norm(step)/trmodel.radius));
                            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                        end
                    end
                end
            end
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
%             interp_error = check_interpolation(trmodel); % TO BE REMOVED!
%             if pred > 0 && rho > eta_2 && pred < interp_error
%                 gamma_dec = gamma_1;
%                 trmodel.radius = gamma_dec*trmodel.radius;
%                 trmodel = improve_model(trmodel, fphi, options);
%             end
            if trmodel.radius < tol_radius
                finish = true;
                break
            end
            history_solution(end+1).x = x;
            history_solution(end).rho = rho;
            history_solution(end).radius = trmodel.radius;
            history_solution(end).px = px;
            history_solution(end).fx = trmodel.fvalues(1,1);
        end
    end
    history_solution(end+1).x = x;
    history_solution(end).rho = rho;
    history_solution(end).radius = trmodel.radius;
    history_solution(end).px = px;
    history_solution(end).fx = trmodel.fvalues(1,1);
    if trmodel.radius < 0.1
        1;
    end
%     try
%         interp_error = check_interpolation(trmodel); % TO BE REMOVED!
%     catch erro
%         retrhow(erro);
%     end
%     if pred > 0 && rho > eta_2 && pred < interp_error
%         gamma_dec = gamma_1;
%         trmodel.radius = gamma_dec*trmodel.radius;
%         trmodel = improve_model(trmodel, fphi, options);
%     end
    if trmodel.radius < tol_radius
        finish = true;
    end
    end

end