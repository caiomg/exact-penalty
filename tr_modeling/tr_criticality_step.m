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

[~, fmodel.g] = get_model_matrices(model, 0);
cmodel = extract_constraints_from_tr_model(model, con_lb, con_ub);
[ind_eactive, ~] = identify_new_constraints(cmodel, epsilon, []);
[Q, R, ind_qr] = update_factorization(cmodel, [], [], ind_eactive, true);
pseudo_gradient = l1_pseudo_gradient(fmodel.g, p_mu, cmodel, ind_qr, true);

% Criticality measure
measure = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, [cmodel(ind_qr).c]');

if norm(measure) <= Lambda
    while true
        [multipliers, tol_multipliers] = ...
            l1_estimate_multipliers(fmodel, cmodel, p_mu, ind_qr, Q, R, x, bl, bu);
        if sum(multipliers < -tol_multipliers | p_mu < multipliers - tol_multipliers)
            [Q, R, N, ind_qr, ind_eactive] = ...
                    l1_drop_constraint(cmodel, Q, R, ind_qr, ind_eactive, p_mu, ...
                                       multipliers, tol_multipliers);
            pseudo_gradient = l1_pseudo_gradient(fmodel.g, p_mu, cmodel, ind_qr, true);
            measure = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, [cmodel(ind_qr).c]');
            if norm(measure) > Lambda
                break
            end
        else
            break
        end
    end
end



while (model.radius > crit_mu*measure)
    model.radius = omega*model.radius;
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

    [~, fmodel.g] = get_model_matrices(model, 0);
    cmodel = extract_constraints_from_tr_model(model, con_lb, con_ub);
    [ind_eactive, ~] = identify_new_constraints(cmodel, ...
                                                      epsilon, []);
    [Q, R, ind_qr] = update_factorization(cmodel, [], [], ...
                                             ind_eactive, true);
    measure = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, ...
                                     bu, [cmodel(ind_qr).c]');
    if norm(measure) <= Lambda
        while true
            [multipliers, tol_multipliers] = ...
                l1_estimate_multipliers(fmodel, cmodel, p_mu, ind_qr, Q, R, x, bl, bu);
            if sum(multipliers < -tol_multipliers | p_mu < multipliers - tol_multipliers)
                [Q, R, N, ind_qr, ind_eactive] = ...
                        l1_drop_constraint(cmodel, Q, R, ind_qr, ind_eactive, p_mu, ...
                                           multipliers, tol_multipliers);
                pseudo_gradient = l1_pseudo_gradient(fmodel.g, p_mu, cmodel, ind_qr, true);
                measure = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, bu, [cmodel(ind_qr).c]');
                if norm(measure) > Lambda
                    break
                end
            else
                break
            end
        end
    end
    pseudo_gradient = l1_pseudo_gradient(fmodel.g, p_mu, cmodel, ind_qr, true);

    measure = l1_criticality_measure(x, pseudo_gradient, Q, R, bl, ...
                                     bu, [cmodel(ind_qr).c]');
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

% The final radius is increased not to make the reduction drastic
model.radius = min(max(model.radius, beta*norm(measure)), initial_radius);
identified = length(ind_eactive);
while true
    [active_larger, ~] = identify_new_constraints(cmodel, epsilon/factor_epsilon, []);
    if identified == length(active_larger) && epsilon/factor_epsilon <= epsilon0
        epsilon = epsilon/factor_epsilon;
    else
        break
    end
end

                       
end