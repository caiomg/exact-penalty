function v = tr_vertical_step_new(fmodel, cmodel, mu, h, ind_eactive, radius, x0, bl, bu)

    if isempty(bl)
        bl = -inf(size(x0));
    end
    if isempty(bu)
        bu = inf(size(x0));
    end
opts_lt.LT = true;
tol_r = 1e-8;

dimension = size(h, 1);
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
    bounds_included = false(dimension, 1);
    Z = zeros(dimension, 1);
    while true
        Qa = orth(An);
        v = -Qa*((A2'*Qa)\uphi); % possibly outside TR
        if norm(v) > r_radius
            v = (v/norm(v))*r_radius;
        end
        xn = x0 + h;
        bl_active = xn <= bl & v < 0;
        bu_active = xn >= bu & v > 0;
        included = 0;
        for k = 1:dimension
           if bl_active(k) && bounds_included(k) == false
               An(k, :) = zeros(1, size(An, 2));
               included = k;
               break % for
           elseif bu_active(k) && bounds_included(k) == false
               An(k, :) = zeros(1, size(An, 2));
               included = k;
               break % for
           end
        end
        if included ~= 0
            bounds_included(included) = true;
        else
            break % while
        end
    end
    pred_hv = predict_descent(h_fmodel, h_constraints, v, mu);
    if true %pred_hv < -0.75*pred_h
        [v, ~, status] = line_search_full_domain(h_fmodel, h_constraints, mu, v, r_radius);
        if ~status
            v = zeros(size(v));
        end
    else
        % Will use as it is
        1;
    end

else
    v = zeros(size(h));
end




end



