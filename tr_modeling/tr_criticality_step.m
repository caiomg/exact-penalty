function [model, epsilon] = tr_criticality_step(model, f, epsilon, mu, ...
                                     options)
% CRITICALITY_STEP -- ensures model is sufficiently poised and with
% a radius comparable to the gradient
%


mu = options.criticality_mu; % factor between radius and
                             % criticality measure
omega = options.criticality_omega; % factor used to reduce radius
beta = options.criticality_beta; % to ensure the final radius
                                 % reduction is not drastic
tol_radius = options.tol_radius; % tolerance of TR algorithm
tol_f = options.tol_f;

initial_radius = model.radius;
model = improve_model(model, f, options);
[~, f_grad] = get_model_matrices(model, 0);
cmodel = extract_constraints_from_tr_model(model);
[ind_eactive, ind_eviolated] = identify_new_constraints(cmodel, epsilon, []);
[N, Q, R, ind_qr] = update_factorization(cmodel, [], [], ind_eactive, true);
pseudo_gradient = l1_pseudo_gradient(f_grad, mu, cmodel, ind_eviolated);
q1 = N'*pseudo_gradient;
q2 = [cmodel(ind_qr).c]';

measure = sqrt(q1'*q1 + q2'*q2);


while (model.radius > mu*measure)
    model.radius = omega*model.radius;
    epsilon = 0.5*epsilon;
    
    model = improve_model(model, f, options);
    [~, f_grad] = get_model_matrices(model, 0);
    cmodel = extract_constraints_from_tr_model(model);
    [ind_eactive, ind_eviolated] = identify_new_constraints(cmodel, ...
                                                      epsilon, []);
    [N, Q, R, ind_qr] = update_factorization(cmodel, [], [], ...
                                             ind_eactive, true);
    pseudo_gradient = l1_pseudo_gradient(f_grad, mu, cmodel, ind_eviolated);
    q1 = N'*pseudo_gradient;
    q2 = [cmodel(ind_qr).c]';

    measure = sqrt(q1'*q1 + q2'*q2);
    if (model.radius < tol_radius || ...
        beta*measure < tol_f)
        % Better break.
        % Not the end of this algorithm, but satisfies stopping
        % condition for outer algorithm anyway...
        break;
    end
end

% The final radius is increased not to make the reduction drastic
model.radius = min(max(model.radius, beta*norm(measure)), initial_radius);
                       
end