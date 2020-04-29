function [x, history_solution] = l1_penalty_solve(f, phi, con_lb, con_ub, initial_points, ...
                                                  mu, epsilon, delta, ...
                                                  Lambda, bl, bu, options)
%L1_PENALTY Summary of this function goes here
%   Detailed explanation goes here

defaultoptions = struct('tol_radius', 1e-6, ...
                        'tol_f', 1e-6, ...
                        'tol_measure', 1e-5, ...
                        'tol_con', 1e-5, ...
                        'eps_c', 1e-4, ...
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
tol_measure = options.tol_measure;
tol_con = options.tol_con;
tol_radius = options.tol_radius;
tol_f = options.tol_f;
initial_radius = options.initial_radius;
epsilon_decrease_radius_threshold = initial_radius;
epsilon_decrease_measure_threshold = 1e3*tol_measure;
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

p = @(x) l1_function(f, phi, con_lb, con_ub, mu, x);



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
history_solution.sigma = nan;
history_solution.ns = 0;
history_solution.mchange = nan;
history_solution.epsilon = epsilon;
history_solution.Lambda = Lambda;
history_solution.pred = nan;
history_solution.ared = nan;
history_solution.polynomial = trmodel.modeling_polynomials{1};
count_inf = 0;
evaluate_step = true;

while ~finish

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
    
    [measure, d, is_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel, cmodel, ...
                                                     x, mu, epsilon, bl, bu);

    if measure > eps_c
        tr_criticality_step_executed = false;
    else
        tr_criticality_step_executed = true;
        [trmodel, epsilon, epsilon_decrease_measure_threshold, ...
         epsilon_decrease_radius_threshold] = tr_criticality_step(trmodel, ...
                                                          fphi, mu, ...
                                                          epsilon, ...
                                                          bl, bu, ...
                                                          con_lb, ...
                                                          con_ub, ...
                                                          epsilon_decrease_measure_threshold, ...
                                                          epsilon_decrease_radius_threshold, ...
                                                          options);

        if ~strcmp(options.basis, 'dummy')
            trmodel.modeling_polynomials = compute_polynomial_models(trmodel);
            [fx, fmodel.g, fmodel.H] = get_model_matrices(trmodel, 0);
            cmodel = extract_constraints_from_tr_model(trmodel, con_lb, con_ub);
        else
            [fx, fmodel.g, fmodel.H] = f(x);
            cmodel = evaluate_constraints(phi, x, con_lb, con_ub);
        end

        [measure, d, is_eactive] = ...
            l1_criticality_measure_and_descent_direction(fmodel, ...
                                                         cmodel, x, ...
                                                         mu, epsilon, bl, bu);
    end

    % Stopping conditions
    if measure < tol_measure
        eactive_norm_inf = norm([cmodel(is_eactive).c], inf);
        if eactive_norm_inf < tol_con
            if max([cmodel.c]) > tol_con
                break % Address constraint violation outside
            else
                break
            end
        end
    end

    geometry_ok = is_lambda_poised(trmodel, options);

    [point_computed, pred, Lambda] = l1_trust_region_step(fmodel, cmodel, x, ...
                                               epsilon, Lambda, mu, ...
                                               trmodel.radius, bl, bu);
    trial_point = project_to_bounds(point_computed, bl, bu);
    if norm(trial_point - point_computed) > 0 
        'Debug this';
    end
    if pred <= 0
        'Debug this';
    end
    s = trial_point - x;
    evaluate_step = pred > 0; % Another logic can be used
    if evaluate_step
        [p_trial, trial_fvalues] = p(trial_point);
        if ~isempty(find(~isfinite(trial_fvalues), 1))
            rho = -inf;
            ared = nan;
        else
            ared = ...
              evaluate_p_descent(trmodel.fvalues(:, trmodel.tr_center), ...
                                 trial_fvalues, con_lb, con_ub, mu);
            rho = ared/pred;
        end
        if rho >= eta_2 || (rho > eta_1 && geometry_ok)
            px = p_trial;
            x = trial_point;
            if strcmp(options.basis, 'dummy')
                mchange_flag = 4;
                trmodel.points_abs(:, 1) = x;
                trmodel.tr_center = 1;
            else
                try
                [trmodel, mchange_flag] = ...
                    change_tr_center(trmodel, trial_point, ...
                                     trial_fvalues, options);
                catch this_error
                    rethrow(this_error);
                end
            end
        elseif isfinite(rho)
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
        ared = 0;
        [trmodel, mchange_flag] = ...
                ensure_improvement(trmodel, fphi, bl, bu, options);
        if strcmp(options.basis, 'dummy')
                mchange_flag = 4;
        end
    end
    if rho < eta_2
        gamma_dec = gamma_1;
        trmodel.radius = gamma_dec*trmodel.radius;
    elseif isfinite(rho)
        s_norm = norm(s, inf);
        s_norm = min(trmodel.radius, s_norm); %small correction
        radius_inc = max(1, gamma_2*(s_norm/trmodel.radius));
        trmodel.radius = min(radius_inc*trmodel.radius, radius_max);
    else
        warning('cmg:infinite_rho', 'Invalid rho = %g', rho);
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
    history_solution(iter).sigma= measure;
    history_solution(iter).ns = norm(s);
    history_solution(iter).mchange = mchange_flag;
    history_solution(iter).epsilon = epsilon;
    history_solution(iter).Lambda= Lambda;
    history_solution(iter).pred = pred;
    history_solution(iter).ared = ared;
    history_solution(iter).criticality_step = tr_criticality_step_executed;

    offending_pivot = find(isinf(trmodel.pivot_values), 1);
    if ~isempty(offending_pivot)
        warning('cmg:offending_pivot', 'Pivot value %g at iter %d, radius: %g\n', ...
                trmodel.pivot_values(offending_pivot), iter, trmodel.radius);
        trmodel = rebuild_model(trmodel, options);
    end

    if trmodel.radius < tol_radius || iter > max_iter
        finish = true;
    end
    if debug_on
        try
            %check_nfp_polynomials(trmodel);
        catch this_error
           rethrow(this_error);
        end
        if iter == inspect_iteration
            warning('cmg:inspect_iteration', 'Iteration %d reached', iter);
        end
    end
end
% In case we decided to preallocate
history_solution = history_solution(1:iter);
