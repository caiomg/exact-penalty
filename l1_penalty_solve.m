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

iter = 1;
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
history_solution.q = nan;
history_solution.ns = 0;
while ~finish

    [~, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
    current_constraints = extract_constraints_from_tr_model(trmodel);

    [ind_eactive, ~] = ...
        identify_new_constraints(current_constraints, epsilon, []);
    [N, Q, R, ind_qr] = update_factorization(current_constraints, ...
                                             Q, R, ind_eactive, false);
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, current_constraints, ...
                                         ind_qr, true);

    q = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, ...
                               [current_constraints(ind_qr).c]');

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
        q = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, ...
                                   [current_constraints(ind_qr).c]');
    end

    if (norm(q) <= max(Lambda, tol_g))
            [multipliers, tol_multipliers] = ...
            l1_estimate_multipliers(fmodel, current_constraints, mu, ind_qr, ...
                                    Q, R, x, bl, bu);

        % Stopping conditions
        if ~sum(multipliers < -tol_multipliers | ...
                mu < multipliers - tol_multipliers) 
            if q < tol_g
                phih = [current_constraints(ind_qr).c]';
                if norm(phih) < tol_con
                    if false %max([current_constraints.c]) > tol_g
                        mu = mu*10
                        % should recalculate px
                    else
                        break
                    end
                end
            end
        end

        % Are there conditions for dropping one constraint?
        if sum(multipliers < -tol_multipliers | mu < multipliers - ...
               tol_multipliers)
            while sum(multipliers < -tol_multipliers | mu < multipliers ...
                      - tol_multipliers)
                [Q, R, N, ind_qr, ind_eactive] = ...
                        l1_drop_constraint(current_constraints, Q, R, ind_qr, ind_eactive, mu, ...
                                           multipliers, tol_multipliers);
                [multipliers, tol_multipliers] = ...
                    l1_estimate_multipliers(fmodel, current_constraints, ...
                                            mu, ind_qr, Q, R, x, bl, bu);
                % break
            end
%             pseudo_gradient = l1_pseudo_gradient(fmodel.g, mu, ...
%                                                  current_constraints, ...
%                                                  ind_qr, true);
            q = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, ...
                                       [current_constraints(ind_qr).c]');
        end
    end

    geometry_ok = is_lambda_poised(trmodel, options);
    if norm(q) > Lambda
        % First-order step
        try
        [h, pred_h] = l1_horizontal_step(fmodel, current_constraints, ...
                                         mu, x, ind_qr, Q, R, ...
                                         trmodel.radius, bl, bu);
        [hc, pred_hc] = l1_horizontal_cauchy_step(fmodel, ...
                                                  current_constraints, ...
                                                  mu, x0, ind_qr, ...
                                                  Q, R, trmodel.radius, ...
                                                  bl, bu);
            if pred_h < pred_hc
                h = hc;
                pred_h = pred_hc;
            end
        catch err1
            rethrow(err1);
        end
        v = tr_vertical_step_new(fmodel, current_constraints, Q, R, ...
                                 mu, h, ind_qr, trmodel.radius, x, bl, bu);

        s = (h + v);
        pred = predict_descent(fmodel, current_constraints, s, mu, []);

        if pred < pred_h
            s = h;
            pred = pred_h;
        end
    else
        % Step including multipliers
        if norm(N'*pseudo_gradient) >= tol_g
            try
                [h, pred_h] = l1_horizontal_step(fmodel, ...
                                                 current_constraints, ...
                                                 mu, x, ind_qr, Q, ...
                                                 R, trmodel.radius, ...
                                                 bl, bu, multipliers);
                [hc, pred_hc] = l1_horizontal_cauchy_step(fmodel, ...
                                                          current_constraints, mu, x, ind_qr, Q, R, trmodel.radius, bl, bu, multipliers);
                if pred_h < pred_hc
                    h = hc;
                end
            catch err1
                rethrow(err1);
            end
            v = tr_vertical_step_new(fmodel, current_constraints, ...
                                     Q, R, mu, h, ind_qr, trmodel.radius, ...
                                     x, bl, bu);
            s = correct_step_to_bounds(x, h + v, bl, bu);
        else
            v = tr_vertical_step_new(fmodel, current_constraints, ...
                                     Q, R, mu, zeros(dimension,1), ind_qr, ...
                                     trmodel.radius, x, bl, bu);
            s = correct_step_to_bounds(x, v, bl, bu);
        end
        pred = predict_descent(fmodel, current_constraints, s, mu, []);
        normphi = norm([current_constraints(ind_eactive).c], 1);
        ppgrad = N'*pseudo_gradient;
    end
    if q < Lambda
        if pred < delta*(norm(ppgrad)^2 + normphi)
            evaluate_step = false;
        else
            evaluate_step = true;
        end
    else
        if pred <= 0 % || (norm(s) < 0.05*trmodel.radius && ~geometry_ok)
            evaluate_step = false;
        else
            evaluate_step = true;
        end
    end
    if evaluate_step
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
            trmodel = move_trust_region(trmodel, x, trial_fvalues, ...
                                   fphi, bl, bu, options);
            px = p_trial;
        else
            trmodel = try_to_add_interpolation_point(trmodel, trial_point, ...
                                                trial_fvalues, ...
                                                fphi, bl, bu, options);
        end
    else
        rho = -inf;
        if ~geometry_ok
            trmodel = improve_model(trmodel, fphi, bl, bu, options);
        end
        if geometry_ok && (q < Lambda || pred <= 0)
            [epsilon, Lambda] = l1_reduce_lambda(epsilon, Lambda, ...
                                                 current_constraints, ...
                                                 fmodel.g, mu, Q, ...
                                                 R, tol_g, tol_con); 
        end
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

    iter = iter + 1;
    history_solution(iter).x = x;
    history_solution(iter).rho = rho;
    history_solution(iter).radius = trmodel.radius;
    history_solution(iter).px = px;
    history_solution(iter).fx = trmodel.fvalues(1,1);
    history_solution(iter).q = q;
    history_solution(iter).ns = norm(s);


    if trmodel.radius < tol_radius
        finish = true;
    end
    if length(history_solution) > 55
        1;
    end
end
% In case we decided to preallocate
history_solution = history_solution(1:iter);

