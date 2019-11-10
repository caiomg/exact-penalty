function [model, epsilon] = tr_criticality_step(model, funcs, p_mu, ...
                                                epsilon, Lambda, bl, ...
                                                bu, con_lb, con_ub, ...
                                                options, one_pass)
% CRITICALITY_STEP -- ensures model is sufficiently poised and with
% a radius comparable to the gradient
%

if nargin < 11 || isempty(one_pass)
    one_pass = false;
end

crit_mu = options.criticality_mu; % factor between radius and
                             % criticality measure
omega = options.criticality_omega; % factor used to reduce radius
beta = options.criticality_beta; % to ensure the final radius
                                 % reduction is not drastic
tol_radius = options.tol_radius; % tolerance of TR algorithm
tol_f = options.tol_f;
tol_con = options.tol_con;
factor_epsilon = 0.5;
epsilon0 = epsilon;

x = model.points_abs(:, model.tr_center);
initial_radius = model.radius;
model_changed = false;
if has_distant_points(model, options) || is_old(model, options)
    model = rebuild_model(model, options);
    model_changed = true;
end
while ~is_lambda_poised(model, options)
    model = ensure_improvement(model, funcs, bl, bu, options);
    model_changed = true;
end
if model_changed
    model.modeling_polynomials = compute_polynomial_models(model);
end

if ~strcmp(options.basis, 'dummy')
    [fx, fmodel.g, fmodel.H] = get_model_matrices(model, 0);
    cmodel = extract_constraints_from_tr_model(model, con_lb, con_ub);
else
    [fx, fmodel.g, fmodel.H] = funcs{1}(x);
    cmodel = evaluate_constraints({funcs{2:end}}, x, con_lb, con_ub);
end

% [ind_eactive, ~] = identify_new_constraints(cmodel, epsilon, []);
% [Q, R, ind_qr] = update_factorization(cmodel, [], [], ind_eactive, false);

% [q1, q2, p_grad, d_drop, Q_drop, R_drop, ind_qr_drop, multipliers] = ...
%    l1_measure_criticality(fmodel, cmodel, mu, Q, R, ind_qr, x, ...
%                           bl, bu, max(Lambda, eps_c));
measure = l1_criticality_measure_and_descent_direction(fmodel, ...
                                                       cmodel, p_mu, ...
                                                       x, epsilon, bl, bu);



while (model.radius > crit_mu*measure)
    model.radius = omega*model.radius;
    %epsilon = max(0.5*tol_con, factor_epsilon*epsilon);
    epsilon = factor_epsilon*epsilon;
    
    model_changed = false;
    if has_distant_points(model, options) || is_old(model, options)
        model = rebuild_model(model, options);
        model_changed = true;
    end
    while ~is_lambda_poised(model, options)
        model = ensure_improvement(model, funcs, bl, bu, options);
        model_changed = true;
    end
    if model_changed
        model.modeling_polynomials = compute_polynomial_models(model);
    end

    if ~strcmp(options.basis, 'dummy')
        [fx, fmodel.g, fmodel.H] = get_model_matrices(model, 0);
        cmodel = extract_constraints_from_tr_model(model, con_lb, con_ub);
    end

    % [ind_eactive, ~] = identify_new_constraints(cmodel, epsilon, []);
    % [Q, R, ind_qr] = update_factorization(cmodel, [], [], ind_eactive, false);

    % [q1, q2, p_grad, d_drop, Q_drop, R_drop, ind_qr_drop, multipliers] = ...
    %    l1_measure_criticality(fmodel, cmodel, mu, Q, R, ind_qr, x, ...
    %                           bl, bu, max(Lambda, eps_c));

    measure = l1_criticality_measure_and_descent_direction(fmodel, ...
                                                           cmodel, p_mu, ...
                                                           x, epsilon, bl, bu);

    
    if (model.radius < tol_radius || ...
        (beta*measure < tol_f && model.radius < 100*tol_radius))
        % Better break.
        % Not the end of this algorithm, but satisfies stopping
        % condition for outer algorithm anyway...
        break
    end
    if one_pass
        break
    end
end

while true

    epsilon_larger = epsilon/factor_epsilon;
    measure_larger = l1_criticality_measure_and_descent_direction(fmodel, ...
                                                       cmodel, p_mu, ...
                                                       x, epsilon_larger, bl, bu);

    if model.radius <= crit_mu*measure_larger...
            && epsilon/factor_epsilon <= epsilon0
        measure = measure_larger;
        epsilon = epsilon_larger;
    else
        break
    end
end
% The final radius is increased not to make the reduction drastic
model.radius = min(max(model.radius, beta*norm(measure)), initial_radius);

                       
end