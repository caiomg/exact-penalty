function v = tr_vertical_step_new(fmodel, cmodel, mu, h, ind_eactive, ind_eviolated, radius)


opts_lt.LT = true;
tol_r = 1e-8;

n_cons = length(cmodel);
ls_steps = 10;
ls_factor = 0.75;
pv = @(s) -predict_descent(fmodel, cmodel, s, mu);
r_radius = sqrt(max(0, radius^2 - h'*h));
n_constraints = length(cmodel);
h_fmodel = shift_model(fmodel, h);
for m = 1:n_constraints
    h_constraints(m) = shift_model(cmodel(m), h);
end

pred_h = predict_descent(fmodel, cmodel, h, mu);

if ~isempty(ind_eactive) && r_radius > tol_r
    [uphi, A2] = update_constraint_information(cmodel, ind_eactive, h);
    
    An = [cmodel(ind_eactive).g]; % Change to include perturbations
    Qa = orth(An);
    v = -Qa*((A2'*Qa)\uphi); % possibly outside TR

    pred_hv = predict_descent(fmodel, cmodel, v, mu);
    if pred_hv < 0.25*pred_h
        [v, ~, status] = line_search_full_domain(h_fmodel, h_constraints, mu, v, r_radius);
        if ~status
            v = zeros(size(v));
        end
    end

else
    v = zeros(size(h));
end




end



