function v = tr_vertical_step_new(fmodel, cmodel, Q, R, mu, h, ind_eactive, radius, x0, bl, bu)

if isempty(bl)
    bl = -inf(size(x0));
end
if isempty(bu)
    bu = inf(size(x0));
end
opts_lt.LT = true;
tol_r = 1e-8;
tol_s = 1e-10;


dimension = size(h, 1);
n_cons = length(cmodel);
ls_steps = 10;
ls_factor = 0.75;

r_columns = size(R, 2);
N = Q(:, r_columns+1:end);

pv = @(s) -predict_descent(fmodel, cmodel, s, mu);
r_radius = sqrt(max(0, radius^2 - h'*h));
n_constraints = length(cmodel);
h_fmodel = shift_model(fmodel, h);
for m = 1:n_constraints
    h_constraints(m) = shift_model(cmodel(m), h);
end

pred_h = predict_descent(fmodel, cmodel, h, mu);

if ~isempty(ind_eactive) && r_radius > tol_r
    hv = h;
    
    An = [cmodel(ind_eactive).g]; % Change to include perturbations
    bounds_included = false(dimension, 1);
    Z = zeros(dimension, 1);
    while true
        [uphi, A2] = update_constraint_information(cmodel, ind_eactive, hv);
        A_ext = [A2, N];
        uphi_ext = [uphi; zeros(size(N, 2), 1)];
        v = -(A_ext'\uphi_ext); % possibly outside TR
        if norm(v) > r_radius
            v = (v/norm(v))*r_radius;
        end
        xn = x0 + hv + v;
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
            lower_breakpoints = (bl - xn)./v;
            upper_breakpoints = (bu - xn)./v;
            bp = min([lower_breakpoints(lower_breakpoints > 0); ...
                      upper_breakpoints(upper_breakpoints > 0)]);
            if isempty(bp) ||norm(bp*v) < tol_s
                break % while
            elseif bp >= 1
                hv = hv + v;
                break
            else
                hv = hv + bp*v;
            end
        end
        break
    end

    v = hv - h;
    while true
        pred_hv = predict_descent(h_fmodel, h_constraints, v, mu);
        if pred_hv < -0.5*abs(pred_h)%(pred_h > 0 && pred_hv < -0.5*pred_h) || ...
                %(pred_h <= 0 && pred_hv < 0)
            v = v/2;
        else
            break
        end
    end

else
    v = zeros(size(h));
end




end



