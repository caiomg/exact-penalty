function model = tr_criticality_step(model, functions, measure, options)
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
model = complete_interpolation_set_mn(model, functions, options);
while (model.radius > mu*measure(model))
    model.radius = omega*model.radius;
    model = update_model_new(model, functions, options);
    model = complete_interpolation_set_mn(model, functions, options);
    if (model.radius < tol_radius || ...
        beta*norm(measure(model)) < tol_f)
        % Better break.
        % Not the end of this algorithm, but satisfies stopping
        % condition for outer algorithm anyway...
        break
    end
end

% The final radius is increased not to make the reduction drastic
model.radius = min(max(model.radius, beta*norm(measure(model))), initial_radius);
                       
end