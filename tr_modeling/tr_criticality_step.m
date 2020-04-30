function [model, epsilon, epsilon_decrease_measure_threshold, ...
          epsilon_decrease_radius_threshold] = tr_criticality_step(model, ...
                                                      funcs, ...
                                                      p_mu, ...
                                                      epsilon, ...
                                                      lb, ...
                                                      ub, ...
                                                      con_lb, ...
                                                      con_ub, ...
                                                      epsilon_decrease_measure_threshold, ...
                                                      epsilon_decrease_radius_threshold, ...
                                                      options)
% CRITICALITY_STEP -- ensures model is sufficiently poised and with
% a radius comparable to the gradient
%

crit_mu = options.criticality_mu; % factor between radius and
                                  % criticality measure
omega = options.criticality_omega; % factor used to reduce radius
beta = options.criticality_beta; % to ensure the final radius
                                 % reduction is not drastic
tol_radius = options.tol_radius; % tolerance of TR algorithm
tol_measure = options.tol_measure;
tol_con = options.tol_con;
factor_epsilon = 0.5;
epsilon0 = epsilon;
beta_3 = 1;

gamma_inc = options.gamma_inc;

x = model.points_abs(:, model.tr_center);
initial_radius = model.radius;
dim = size(x, 1);
model_changed = false;
if has_distant_points(model, options) || is_old(model, options)
    model = rebuild_model(model, options);
    model_changed = true;
end
change_count = 0;
while ~is_lambda_poised(model, options)
    [model, mchange_flag] = ensure_improvement(model, funcs, lb, ub, options);
    if mchange_flag == 4
        change_count = change_count + 1;
        if change_count > dim
            model.radius = omega*model.radius;
            change_count = 0;
            if has_distant_points(model, options) || is_old(model, options)
                model = rebuild_model(model, options);
            end
        end
    end
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

[measure, ~, is_eactive] = ...
    l1_criticality_measure_and_descent_direction(fmodel, cmodel, x, ...
                                                 p_mu, epsilon, lb, ub);
eactive_norm = norm([cmodel(is_eactive).c], inf);

detected_convergence_of_main_algorithm = false;
while (model.radius > crit_mu*measure)

    if measure < epsilon_decrease_measure_threshold ...
            && model.radius < epsilon_decrease_radius_threshold ...
            && eactive_norm > beta_3*measure
        epsilon = factor_epsilon*epsilon;
        epsilon_decrease_measure_threshold = 0.5*epsilon_decrease_measure_threshold;
        epsilon_decrease_radius_threshold = omega*epsilon_decrease_radius_threshold;
    else
        model.radius = omega*model.radius;
    end
    
    model_changed = false;
    if has_distant_points(model, options) || is_old(model, options)
        model = rebuild_model(model, options);
        model_changed = true;
    end
    while ~is_lambda_poised(model, options)
        model = ensure_improvement(model, funcs, lb, ub, options);
        model_changed = true;
    end
    if model_changed
        model.modeling_polynomials = compute_polynomial_models(model);
    end

    if ~strcmp(options.basis, 'dummy')
        [fx, fmodel.g, fmodel.H] = get_model_matrices(model, 0);
        cmodel = extract_constraints_from_tr_model(model, con_lb, con_ub);
    end

    [measure, ~, is_eactive] = ...
        l1_criticality_measure_and_descent_direction(fmodel, cmodel, ...
                                                     x, p_mu, epsilon, lb, ub);
    eactive_norm = norm([cmodel(is_eactive).c], inf);

    
    if (model.radius*gamma_inc < tol_radius || ...
        (measure < tol_measure && eactive_norm < tol_con ...
         && model.radius < 100*tol_radius))
        % Better break.
        % Not the end of this algorithm, but satisfies stopping
        % condition for outer algorithm anyway...
        detected_convergence_of_main_algorithm = true;
        break
    end
end

if ~detected_convergence_of_main_algorithm

    % Epsilon is increased back if it does not spoil the criterion
    % of criticality step
    while true

        epsilon_larger = epsilon/factor_epsilon;
        if epsilon_larger <= epsilon0
            
            measure_larger = l1_criticality_measure_and_descent_direction(fmodel, ...
                                                          cmodel, x, p_mu, ...
                                                          epsilon_larger, lb, ub);
            if model.radius <= crit_mu*measure_larger
                measure = measure_larger;
                epsilon = epsilon_larger;
                epsilon_decrease_measure_threshold = 2*epsilon_decrease_measure_threshold;
                epsilon_decrease_radius_threshold = epsilon_decrease_radius_threshold/omega;
            else
                break
            end
        else
            break
        end
    end
    % The final radius is increased not to make the reduction drastic
    model.radius = min(max(model.radius, beta*norm(measure)), initial_radius);
end
                       
end