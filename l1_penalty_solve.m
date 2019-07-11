function [x, history_solution] = l1_penalty_solve(f, phi, con_lb, con_ub, initial_points, ...
                                                  mu, epsilon, delta, ...
                                                  Lambda, bl, bu, options)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

defaultoptions = struct('tol_radius', 1e-6, 'tol_f', 1e-6, ...
                        'tol_con', 1e-5, ...
                        'eps_c', 1e-5, ...
                        'eta_1', 0, ...
                        'eta_2', 0.1, ...
                        'pivot_threshold', 1/16, ...
                        'add_threshold', 100, ...
                        'exchange_threshold', 1000,  ...
                        'initial_radius', 1, ...
                        'radius_max', 1e3, ...
                        'radius_factor', 6, ...
                        'radius_factor_extra_tol', 2, ...
                        'gamma_inc', 2, ...
                        'gamma_dec', 0.5, ...
                        'criticality_mu', 50, ...
                        'criticality_beta', 10, ...
                        'criticality_omega', 0.5, ...
                        'max_iter', 1e6, ...
                        'divergence_threshold', 1e10, ...
                        'basis', 'unused option', ...
                        'debug', false, ...
                        'inspect_iteration', 10);

option_names = fieldnames(defaultoptions);
for k = 1:length(option_names)
    if ~isfield(options, option_names(k))
        options.(option_names{k}) = defaultoptions.(option_names{k});
    end
end
n_functions = 1 + numel(phi);
if ~isfield(options, 'model_type')
    for k = n_functions:-1:1
        model_type{k} = 'minimum norm hessian';
    end
    options.model_type = model_type;
end
assert (n_functions == numel(options.model_type));

debug_on = options.debug;
if debug_on
    inspect_iteration = options.inspect_iteration;
else
    inspect_iteration = inf;
end
                    
gamma_1 = options.gamma_dec;
gamma_2 = options.gamma_inc;
eps_c = options.eps_c;
ut_option.UT = true;
eta_1 = options.eta_1;
eta_2 = options.eta_2;
eps_c = options.eps_c;
initial_radius = options.initial_radius;
rel_pivot_threshold = options.pivot_threshold;
max_iter = options.max_iter;

fphi = {f, phi{:}}';

[dim, n_initial_points] = size(initial_points);
if  (~isempty(bl) && ~isempty(find(initial_points(:, 1) < bl, 1))) || ...
        (~isempty(bu) && ~isempty(find(initial_points(:, 1) > bu, 1)))
     % Replace
     initial_points(:, 1) = project_to_bounds(initial_points(:, 1), bl, bu);
end
if n_initial_points == 1
    % Finding a random second point
    old_seed = rng('default');
    second_point = rand(size(initial_points));
    rng(old_seed);
    while norm(second_point, inf) < rel_pivot_threshold
        % Second point must not be too close
        second_point = 2*second_point;
    end
    second_point = (second_point - 0.5)*initial_radius;
    second_point = initial_points(:, 1) + second_point;
    if ~isempty(bu)
        second_point = min(bu, second_point);
    end
    if ~isempty(bl)
        second_point = max(bl, second_point);
    end
    initial_points(:, 2) = second_point;
    n_initial_points = 2;
end
initial_fvalues = zeros(n_functions, n_initial_points);
% Calculating function values for other points of the set
for nf = 1:n_functions
    for k = 1:n_initial_points
        initial_points(:, k) = project_to_bounds(initial_points(:, k), bl, bu);
        initial_fvalues(nf, k) = fphi{nf}(initial_points(:, k));
    end
end


% Initializing model structure
trmodel = tr_model(initial_points, initial_fvalues, initial_radius);
trmodel = rebuild_model(trmodel, options);
trmodel.modeling_polynomials = compute_polynomial_models(trmodel);




linsolve_opts.UT = true;

x0 = initial_points(:, 1);
x = x0;
dim = size(x, 1);
tol_g = 1e-5;
tol_con = options.tol_con;
tol_radius = options.tol_radius;
tol_f = options.tol_f;

p = @(x) l1_function(f, phi, con_lb, con_ub, mu, x);
% QR decomposition of constraints gradients matrix A
Q = zeros(dim, 0);
R = zeros(0, 0);


radius_max = options.radius_max;

exchange_counts = 0;
iter = 1;
finish = false;

if ~strcmp(options.basis, 'dummy')
    [fx, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
    cmodel = extract_constraints_from_tr_model(trmodel, ...
                                                      con_lb, con_ub);
else
    [fx, fmodel.g, fmodel.H] = f(x);
    cmodel = evaluate_constraints(phi, x, con_lb, con_ub);
end
px = fx + mu*sum(max(0, [cmodel.c]));
history_solution.x = x;
rho = nan;
history_solution.rho = rho;
history_solution.radius = trmodel.radius;
history_solution.px = px;
history_solution.fx = fx;
history_solution.q1 = nan;
history_solution.q2 = nan;
history_solution.q3 = nan;
history_solution.ns = 0;
history_solution.mchange = nan;
history_solution.step_type = nan;
history_solution.epsilon = epsilon;
history_solution.Lambda = Lambda;
history_solution.pred = nan;
history_solution.polynomial = trmodel.modeling_polynomials{1};
count_inf = 0;
evaluate_step = true;

global solucao

while ~finish

%     if debug_on
%         check_interpolation(trmodel);
%     end

    if abs(px) > options.divergence_threshold
        break
    end

    if ~strcmp(options.basis, 'dummy')
        trmodel.modeling_polynomials = compute_polynomial_models(trmodel);
        [fx, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
        cmodel = extract_constraints_from_tr_model(trmodel, con_lb, con_ub);
    else
        [fx, fmodel.g, fmodel.H] = f(x);
        cmodel = evaluate_constraints(phi, x, con_lb, con_ub);
    end
    history_solution(iter).polynomial = trmodel.modeling_polynomials{1};
    if iter > numel(solucao) ...
            || norm(history_solution(iter).polynomial.coefficients ...
            - solucao(iter).polynomial.coefficients, inf) ~= 0
        1;
    end

    [q1, q2, q3, d, Q, R, ind_qr, multipliers, lb_active, ub_active] = ...
        l1_measure_criticality_new(fmodel, cmodel, mu, epsilon, x, ...
                                   bl, bu, max(Lambda, eps_c));
%     [q1_check, q2_check, q3_check, d_check, Q_check, R_check, ind_qr_check, multipliers_check, lb_active_check, ub_active_check] = ...
%         l1_measure_criticality_new(fmodel, cmodel, mu, epsilon, x, ...
%                                    bl, bu, max(Lambda, eps_c));
% 	if q1 - q1_check ~= 0 ...
%             || q2 - q2_check ~= 0 ...
%             || q3 - q3_check ~= 0 ...
%             || norm(d - d_check, inf) ~= 0 ...
%             || norm(Q - Q_check, inf) ~= 0 ...
%             || norm(R - R_check, inf) ~= 0
%         warning('cmg:check-difference', 'fumou');
%     end
    if numel(ind_qr) ~= size(R, 2)
        1;
    end

    if max([q1, q2, q3]) > max(eps_c, Lambda)
        tr_criticality_step_executed = false;
    else
        tr_criticality_step_executed = true;
        [trmodel, epsilon] = tr_criticality_step(trmodel, fphi, mu, ...
                                                 epsilon, Lambda, ...
                                                 bl, bu, con_lb, ...
                                                 con_ub, options);

        if ~strcmp(options.basis, 'dummy')
            [fx, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
            cmodel = ...
                extract_constraints_from_tr_model(trmodel, con_lb, con_ub);
        end

        [q1, q2, q3, d, Q, R, ind_qr, multipliers, lb_active, ub_active] = ...
            l1_measure_criticality_new(fmodel, cmodel, mu, epsilon, x, ...
                                       bl, bu, max(Lambda, eps_c));
        while Lambda > max([q1, q2, q3]) && (q2 ~= 0 || q3 ~= 0)
            Lambda = 0.5*Lambda;
        end
    end
    
    if numel(ind_qr) ~= size(R, 2)
        1;
    end

    % Stopping conditions
    if q1 < tol_g && q2 == 0 && q3 == 0
        phih = [cmodel(ind_qr).c]';
        if norm(phih) < tol_con
            if max([cmodel.c]) > tol_con && false
                break % Address this outside
                mu = mu*10;
                p = @(x) l1_function(f, phi, mu, x);
                px = trmodel.fvalues(1, 1) + mu*sum(max(0, trmodel.fvalues(2:end, 1)));
            else
                break
            end
        end
    end

    geometry_ok = is_lambda_poised(trmodel, options);

    if q1 > Lambda && q3 == 0
        % Regular step, away from stationary point
        step_type = 1;
        [s, pred] = iterated_steps_new(fmodel, cmodel, trmodel.radius, ...
                                         mu, x, Q, R, ind_qr, bl, ...
                                         bu, lb_active, ub_active);
        if pred <= 0
            warning('cmg:fumou', 'check this');
        end
    elseif q2 > Lambda || q3 > 0
        % Drop step
        if q2 > q3
            step_type = 2;
        else
            step_type = 3;
        end
        h = null_space_conjugate_gradient(fmodel, cmodel, mu, Q, R, ...
                                          ind_qr, x, d, trmodel.radius, ...
                                          bl, bu, zeros(dim, 1), ...
                                          lb_active, ub_active, false);
        v = l1_range_step(fmodel, cmodel, Q, R, mu, h, ind_qr, ...
                          trmodel.radius, x, bl, bu);
        s = h + v;
        pred = predict_descent(fmodel, cmodel, s, mu, []);
        if pred <= 0
            warning('cmg:fumou', 'check this');
        end

    elseif q2 == 0 % All multipliers in correct range
        % Near stationary point
        % Step including multipliers
        step_type = 4;
        if q1 > 0
            h = null_space_step_complete(fmodel, cmodel, mu, x, ...
                                         ind_qr, Q, R, trmodel.radius, ...
                                         bl, bu, multipliers);

            v = l1_range_step(fmodel, cmodel, Q, R, mu, h, ind_qr, ...
                             trmodel.radius, x, bl, bu);
            s = h + v;
        else
            h = zeros(dim, 1);
            v = l1_range_step(fmodel, cmodel, Q, R, mu, h, ind_qr, ...
                             trmodel.radius, x, bl, bu);
            s = v;
        end
        pred = predict_descent(fmodel, cmodel, s, mu, []);
        normphi = norm([cmodel(ind_qr).c], 1);
        if pred < delta*(q1^2 + normphi)
            Lambda = 0.5*Lambda;
            % Regular step, away from stationary point
            step_type = 5;
            [s, pred] = iterated_steps_new(fmodel, cmodel, trmodel.radius, ...
                                           mu, x, Q, R, ind_qr, bl, ...
                                           bu, lb_active, ub_active);
            if pred <= 0
                warning('cmg:fumou', 'check this');
            end
        end
    else
        % This shouldn't happen
        h = zeros(dim, 1);
        v = zeros(dim, 1);
        s = zeros(dim, 1);
        pred = 0;
        error('cmg:runtime_error', 'Low criticality. Debug this');
    end
    if (pred < tol_f*abs(px)*1e-6 ...
            && norm(s) < min(tol_radius, 0.1*trmodel.radius))% ...
         %|| (pred < tol_f*abs(px)*1e-6) % ...
        % || (pred < tol_radius*1e-2) ...
        evaluate_step = false;
    else
        evaluate_step = true;
    end
    if norm(s) - 1.1*trmodel.radius > 0
%        warning('long step'); 
    end
    if norm(s) < 0.5*tol_radius
%        warning('short step'); 
    end
    if evaluate_step
        trial_point = project_to_bounds(x + s, bl, bu);
        [p_trial, trial_fvalues] = p(trial_point);
        if find(~isfinite(trial_fvalues), 1)
            rho = -inf;
        else
            ared = ...
              evaluate_p_descent(trmodel.fvalues(:, trmodel.tr_center), ...
                                 trial_fvalues, con_lb, con_ub, mu);
            rho = ared/pred;
        end
        if rho > eta_2 || (rho > eta_1 && geometry_ok)
            x = trial_point;
            if strcmp(options.basis, 'dummy')
                mchange_flag = 4;
                trmodel.points_abs(:, 1) = x;
                trmodel.tr_center = 1;
            else
                [trmodel, mchange_flag] = ...
                    change_tr_center(trmodel, trial_point, ...
                                     trial_fvalues, options);
            end
            px = p_trial;
        elseif rho ~= -inf
            if strcmp(options.basis, 'dummy')
                mchange_flag = 4;
            else
                [trmodel, mchange_flag] = ...
                    try_to_add_point(trmodel, trial_point, ...
                                     trial_fvalues, fphi, bl, bu, options);
            end
        else
            if strcmp(options.basis, 'dummy')
                mchange_flag = 4;
            else
                [trmodel, mchange_flag] = ...
                    ensure_improvement(trmodel, fphi, bl, bu, options);
            end
        end
    else % evaluate step
        rho = -inf;
        [trmodel, mchange_flag] = ...
                ensure_improvement(trmodel, fphi, bl, bu, options);
        if strcmp(options.basis, 'dummy')
                mchange_flag = 4;
        end
    end
    if rho <= eta_2 && (mchange_flag == 4 || tr_criticality_step_executed)
        gamma_dec = gamma_1;
        trmodel.radius = gamma_dec*trmodel.radius;
    else
        if rho > eta_2
            radius_inc = max(1, gamma_2*(norm(s)/trmodel.radius));
            trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
        end
    end
    if mchange_flag == 2 && rho >= eta_2
        exchange_counts = exchange_counts + 1;
        if exchange_counts > 2*dim
            [trmodel, mchange_flag] = ensure_improvement(trmodel, fphi, bl, bu, options);
            exchange_counts = 0;
        end
    else
        exchange_counts = 0;
    end

    iter = iter + 1;
    history_solution(iter).x = x;
    history_solution(iter).rho = rho;
    history_solution(iter).radius = trmodel.radius;
    history_solution(iter).px = px;
    history_solution(iter).fx = trmodel.fvalues(1, trmodel.tr_center);
    history_solution(iter).q1 = q1;
    history_solution(iter).q2 = q2;
    history_solution(iter).q3 = q3;
    history_solution(iter).ns = norm(s);
    history_solution(iter).mchange = mchange_flag;
    history_solution(iter).step_type= step_type;
    history_solution(iter).epsilon = epsilon;
    history_solution(iter).Lambda= Lambda;
    history_solution(iter).pred = pred;


    if trmodel.radius < tol_radius || iter > max_iter
        finish = true;
    end
    if abs(trmodel.pivot_values(find(abs(trmodel.pivot_values) > 0, 1, 'last'))) < gamma_1*options.pivot_threshold*min(1, trmodel.radius)
        1;
    end
    if pred < 0 && ~isinf(rho) || norm(s) - 1.1*history_solution(iter-1).radius > 0
        1;
    end
    if debug_on
        try
            check_nfp_polynomials(trmodel);
        catch this_error
           1;%rethrow(this_error);
        end
        if iter == inspect_iteration
            warning('cmg:inspect_iteration', 'Iteration %d reached', iter);
        end
    end
end
% In case we decided to preallocate
history_solution = history_solution(1:iter);

