function [x, history_solution] = l1_penalty_solve(f, phi, initial_points, ...
                                                  mu, epsilon, delta, ...
                                                  Lambda, bl, bu, options)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

defaultoptions = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                       'eps_c', 1e-5, 'eta_1', 0, 'eta_2', 0.1, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 0.5, 'radius_max', 1e3, ...
                        'criticality_mu', 50, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'full quadratic', ...
                        'pivot_threshold', 0.1, 'poised_radius_factor', 2, ...
                        'pivot_imp', 1.1);

option_names = fieldnames(defaultoptions);
for k = 1:length(option_names)
    if ~isfield(options, option_names(k))
        options.(option_names{k}) = defaultoptions.(option_names{k});
    end
end
                    
                    
gamma_0 = 0.0625;
gamma_1 = 0.5;
gamma_2 = 2;
ut_option.UT = true;
eta_1 = options.eta_1;
eta_2 = options.eta_2;

all_f = {f, phi{:}};
n_functions = size(all_f, 2);
[dimension, n_initial_points] = size(initial_points);
initial_fvalues = zeros(n_functions, n_initial_points);
% Calculating function values for other points of the set
for nf = 1:n_functions
    for k = 1:n_initial_points
        initial_points(:, k) = project_to_bounds(initial_points(:, k), bl, bu);
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
  case 'dummy'
    basis = [];
end

% Initializing model structure
trmodel.points = initial_points;
trmodel.fvalues = initial_fvalues;
trmodel.radius = options.initial_radius;
trmodel.basis = basis;

fphi = {f, phi{:}}';

% Completing set of interpolation points and calculating polynomial
% model
while true
    % Recomplete interpolation set and calculate new model
    [trmodel, exitflag] = complete_interpolation_set(trmodel, fphi, bl, bu, options);
    if exitflag >= 0 || trmodel.radius < options.tol_radius
        break
    else
        trmodel.radius = 0.5*trmodel.radius;
    end
end

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

    [ind_eactive, ~] = ...
        identify_new_constraints(current_constraints, epsilon, []);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                             Q, R, ind_eactive, false);
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                         ind_qr, true);

    [~, q1] = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, [current_constraints(ind_qr).c]');
    if (norm(q1) > max(Lambda, tol_g))
        x_prev = x;


        geometry_ok = is_lambda_poised(trmodel, options);
        
        try
        [h1, pred_h1] = l1_horizontal_step(fmodel, current_constraints, mu, x, ind_qr, Q, R, trmodel.radius, bl, bu);
        [hc, pred_hc] = l1_horizontal_cauchy_step(fmodel, current_constraints, mu, ...
                                                  x0, ind_eactive, Q, ...
                                                  R, trmodel.radius, bl, bu);
            if pred_h1 < pred_hc
                h1 = hc;
            end
        catch err1
            rethrow(err1);
        end
        v1 = tr_vertical_step_new(fmodel, current_constraints, Q, R, mu, ...
                                  h1, ind_qr, trmodel.radius, x, bl, bu);

        s = (h1 + v1);
        
%         s = line_search_full_domain(fmodel, current_constraints, mu, s, trmodel.radius);

        pred = predict_descent(fmodel, current_constraints, s, mu, []);
        if pred <= 0 || (norm(s) < 0.0625*trmodel.radius && ~geometry_ok)
            rho = -inf;
            if geometry_ok
                trmodel.radius = 0.5*trmodel.radius;
            else
                trmodel = improve_model(trmodel, fphi, bl, bu, options);
            end
        else
            trial_point = project_to_bounds(x + s, bl, bu);
            [p_trial, trial_fvalues] = p(trial_point);
            if find(isnan(trial_fvalues), 1)
                rho = -inf;
            else
                ared = px - p_trial; % There are better ways
                rho = ared/pred;
            end

            if rho > eta_2 || (rho > eta_1 && geometry_ok)
                x = trial_point;
                step_accepted = true;
                trmodel = move_trust_region(trmodel, x, trial_fvalues, ...
                                       fphi, bl, bu, options);
                px = p_trial;
            else
                step_accepted = false;
                trmodel = try_to_add_interpolation_point(trmodel, trial_point, ...
                                                    trial_fvalues, ...
                                                    fphi, bl, bu, options);
            end
            if rho < eta_2 && geometry_ok
                gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                trmodel.radius = gamma_dec*trmodel.radius;
            else
                if rho > eta_2
                    radius_inc = max(1, gamma_2*(norm(s)/trmodel.radius));
                    trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                end
            end
        end
    else
        while true
          	x_prev = x;
            [~, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
            current_constraints = extract_constraints_from_tr_model(trmodel);
            pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, ...
                                                 current_constraints, ...
                                                 ind_qr, true);

            [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                                  Q, R, ind_eactive, true);

                q = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, [current_constraints(ind_qr).c]');
                if (q < 100*tol_g)
                    [trmodel, epsilon] = tr_criticality_step(trmodel, fphi, epsilon, ...
                                                  mu, bl, bu, options);
                    [~, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
                    current_constraints = ...
                        extract_constraints_from_tr_model(trmodel);
                    [ind_eactive, ~] = ...
                        identify_new_constraints(current_constraints, epsilon, []);
                    [N, Q, R, ind_qr] = update_factorization(current_constraints, [], ...
                                                             [], ...
                                                             ind_eactive, true);
                    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, ...
                                                         current_constraints, ...
                                                         ind_qr, true);
                    q = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, [current_constraints(ind_qr).c]');
                    [multipliers, tol_multipliers] = l1_estimate_multipliers(fmodel, current_constraints, mu, ind_qr, Q, R, x, bl, bu);

                    if ~sum(multipliers < -tol_multipliers |...
                            mu < multipliers - tol_multipliers) 
                        phih = [current_constraints(ind_qr).c]';
                        if q < tol_g
                            if norm(phih) < tol_con
                                if false %max([current_constraints.c]) > tol_g
                                    mu = mu*10
                                    % should recalculate px
                                else
                                    finish = true;
                                    break
                                end
                            end
                        end
                    end
                end
            % calculate multipliers
            [multipliers, tol_multipliers] = l1_estimate_multipliers(fmodel, current_constraints, mu, ind_qr, Q, R, x, bl, bu);

            dropping_constraint = false;
            % Are there conditions for dropping one constraint?
            if sum(multipliers < -tol_multipliers | mu < multipliers - tol_multipliers)
                dropping_constraint = true;
                while sum(multipliers < -tol_multipliers | mu < multipliers - tol_multipliers)
                    [Q, R, N, ind_qr] = ...
                            l1_drop_constraint(Q, R, N, ind_qr, mu, ...
                                               multipliers, tol_multipliers);
                    [multipliers, tol_multipliers] = l1_estimate_multipliers(fmodel, current_constraints, mu, ind_qr, Q, R, x, bl, bu);
                    % break
                end

                geometry_ok = is_lambda_poised(trmodel, options);
                [h, pred_h] = l1_horizontal_step(fmodel, current_constraints, mu, x, ind_qr, Q, R, trmodel.radius, bl, bu, multipliers);
                [hc, pred_hc] = l1_horizontal_cauchy_step(fmodel, current_constraints, mu, x, ind_qr, Q, R, trmodel.radius, bl, bu, multipliers);
                if pred_h < pred_hc
                    h = hc;
                end
                v = tr_vertical_step_new(fmodel, current_constraints, Q, R, mu, h, ind_qr, trmodel.radius, x, bl, bu);
                % h = line_search_full_domain(fmodel, current_constraints, mu, h, trmodel.radius);

                % s = project_to_bounds(x + h, bl, bu) - x;
                s = correct_step_to_bounds(x, h + v, bl, bu);
                pred = predict_descent(fmodel, current_constraints, s, mu, []);
                if pred > delta
                    dropping_succeeded = true;
                    trial_point = project_to_bounds(x + s, bl, bu);
                    [p_trial, trial_fvalues] = p(trial_point);
                    if find(isnan(trial_fvalues), 1)
                        rho = -inf;
                    else
                        ared = px - p_trial; % There are better ways
                        rho = ared/pred;
                    end
                    
                    if rho > eta_2 || (rho > eta_1 && geometry_ok)
                        x = trial_point;
                        step_accepted = true;
                        px = p_trial;
                        trmodel = move_trust_region(trmodel, x, ...
                                                    trial_fvalues, ...
                                                    fphi, bl, bu, options);
                    else
                        step_accepted = false;
                        trmodel = try_to_add_interpolation_point(trmodel, ...
                                                                 trial_point, trial_fvalues, fphi, bl, bu, options);
                    end
                    if rho < eta_2 && geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        if rho > eta_2
                            radius_inc = max(1, gamma_2*(norm(s)/trmodel.radius));
                            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                        end
                    end
                else % min decrease not satisfied
                    dropping_succeeded = false;
                    step_accepted = false;
                    if geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(s)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        trmodel = improve_model(trmodel, fphi, bl, bu, options);
                    end
                end
            else

                geometry_ok = is_lambda_poised(trmodel, options);
                if norm(N'*pseudo_gradient) >= tol_g
%                     d1 = -model.g*(trmodel.radius/norm(model.g));
                    try
                        [h1, pred_h] = l1_horizontal_step(fmodel, current_constraints, mu, x, ind_qr, Q, R, trmodel.radius, bl, bu, multipliers);
                        [hc, pred_hc] = l1_horizontal_cauchy_step(fmodel, current_constraints, mu, x, ind_qr, Q, R, trmodel.radius, bl, bu, multipliers);
                        if pred_h < pred_hc
                            h1 = hc;
                        end
                    catch err1
                        rethrow(err1);
                    end
                    v1 = tr_vertical_step_new(fmodel, current_constraints, Q, R, mu, h1, ind_qr, trmodel.radius, x, bl, bu);
                    % s = project_to_bounds(x + (h1 + v1), bl, bu) - x;
                    s = correct_step_to_bounds(x, h1 + v1, bl, bu);
                else
                    v = tr_vertical_step_new(fmodel, current_constraints, ...
                                             Q, R, mu, zeros(dimension,1), ind_qr, ...
                                             trmodel.radius, x, bl, bu);
                    %s = project_to_bounds(x + v, bl, bu) - x;
                    s = correct_step_to_bounds(x, v, bl, bu);
                end
                pred = predict_descent(fmodel, current_constraints, s, mu, []);
                normphi = norm([current_constraints(ind_eactive).c], 1);
                ppgrad = N'*pseudo_gradient;
                if pred < delta*(norm(ppgrad)^2 + normphi)
                    % Better not to try now
                    rho = -inf;
                    step_accepted = false;
                    if geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(Ns+v)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        trmodel = improve_model(trmodel, fphi, bl, bu, options);
                    end
                else
                    % Compute ared and all...
                    trial_point = project_to_bounds(x + s, bl, bu);
                    [p_trial, trial_fvalues] = p(trial_point);
                    if find(isnan(trial_fvalues), 1)
                        rho = -inf;
                    else
                        ared = px - p_trial; % There are better ways
                        rho = ared/pred;
                    end
                    
                    
                    
                    if rho > eta_2 || (rho > eta_1 && geometry_ok)
                        x = trial_point;
                        step_accepted = true;
                        px = p_trial;
                        trmodel = move_trust_region(trmodel, x, ...
                                                    trial_fvalues, fphi, bl, bu, options);
                    else
                        step_accepted = false;
                        trmodel = try_to_add_interpolation_point(trmodel, trial_point, ...
                                                                 trial_fvalues, ...
                                                                 fphi, bl, bu, options);
                    end
                    if rho < eta_2 && geometry_ok
                        gamma_dec = gamma_1;%max(gamma_0, gamma_1*norm(step)/trmodel.radius);
                        trmodel.radius = gamma_dec*trmodel.radius;
                    else
                        if rho > eta_2
                            radius_inc = max(1, gamma_2*(norm(s)/trmodel.radius));
                            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
                        end
                    end
                end
            end
            if geometry_ok && ...
                    (~step_accepted || ...
                     (dropping_constraint && ~dropping_succeeded))
                [epsilon, Lambda, N, Q, R, ind_eactive, ~, ...
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


    if trmodel.radius < tol_radius
        finish = true;
    end
    if length(history_solution) > 1000
        1;
    end
end

end