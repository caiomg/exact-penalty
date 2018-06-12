function [model, epsilon] = tr_criticality_step(model, f, epsilon, p_mu, ...
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
tol_f = options.tol_f;
factor_epsilon = 0.5;
epsilon0 = epsilon;

initial_radius = model.radius;
model = improve_model(model, f, options);
[~, f_grad] = get_model_matrices(model, 0);
cmodel = extract_constraints_from_tr_model(model);
epsilon = factor_epsilon*epsilon;
[ind_eactive, ind_eviolated] = identify_new_constraints(cmodel, epsilon, []);
[N, Q, R, ind_qr] = update_factorization(cmodel, [], [], ind_eactive, true);

pseudo_gradient = l1_pseudo_gradient(f_grad, p_mu, cmodel, ind_eviolated);
q1 = N'*pseudo_gradient;
q2 = [cmodel(ind_qr).c]';

measure = sqrt(q1'*q1 + q2'*q2);


while (model.radius > crit_mu*measure)
    model.radius = omega*model.radius;
    epsilon = factor_epsilon*epsilon;
    
    model = improve_model(model, f, options);
    [~, f_grad] = get_model_matrices(model, 0);
    cmodel = extract_constraints_from_tr_model(model);
    [ind_eactive, ind_eviolated] = identify_new_constraints(cmodel, ...
                                                      epsilon, []);
    [N, Q, R, ind_qr] = update_factorization(cmodel, [], [], ...
                                             ind_eactive, true);
    pseudo_gradient = l1_pseudo_gradient(f_grad, p_mu, cmodel, ind_eviolated);
    q1 = N'*pseudo_gradient;
    if isempty(q1)
        q1 = 0;
    end
    q2 = [cmodel(ind_qr).c]';
    if isempty(q2)
        q2 = 0;
    end

    measure = sqrt(q1'*q1 + q2'*q2);
    try
    if (model.radius < tol_radius || ...
        (beta*measure < tol_f && model.radius < 100*tol_radius))
        % Better break.
        % Not the end of this algorithm, but satisfies stopping
        % condition for outer algorithm anyway...
        break;
    end
    catch erro
        rethrow(erro)
    end
end

% The final radius is increased not to make the reduction drastic
model.radius = min(max(model.radius, beta*norm(measure)), initial_radius);
identified = length(ind_eactive);
while true
    [ind_eactive2, ~] = identify_new_constraints(cmodel, epsilon/factor_epsilon, []);
    if identified == length(ind_eactive2) && epsilon/factor_epsilon <= epsilon0
        epsilon = epsilon/factor_epsilon;
    else
        break
    end
end

                       
end