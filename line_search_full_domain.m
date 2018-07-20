function [s, pred, status] = line_search_full_domain(fmodel, cmodel, mu, s0, radius)

if nargin < 5 || isempty(radius)
    radius = norm(s);
end
dim = size(fmodel.g, 1);

pred0 = predict_descent(fmodel, cmodel, s0, mu, []);
if pred0 < 0
    n_constraints = length(cmodel);
    active = [];
    gv = zeros(dim, 1);
    Bv = zeros(dim);
    for n = 1:n_constraints
        if cmodel(n).c > 0 || ...
               (cmodel(n).c == 0 && ...
                (cmodel(n).g'*s0 > 0 || ...
                 (cmodel(n).g'*s0 == 0 && s0'*cmodel(n).H*s0 > 0)))
            gv = gv + cmodel(n).g;
            Bv = Bv + cmodel(n).H;
            active(end+1) = n;
        end
    end
    g = fmodel.g + mu*gv;
    B = fmodel.H + mu*Bv;
    alpha_s = 0;
    gamma = [];
    total_descent = 0;
    IM = [];

    % Calculating TR constraint
    if isinf(radius)
        gamma_tr = inf;
    else
        bpoint = roots([s0'*s0, 0, - radius^2]);
        gamma_tr = min(bpoint(bpoint > 0));
        if isempty(gamma_tr) || ~isreal(gamma_tr)
            error;
        end
    end
    for n = 1:n_constraints
        Hc = cmodel(n).H;
        gc = cmodel(n).g;
        cc = cmodel(n).c;
        bpoint = roots([0.5*(s0'*Hc*s0),  gc'*s0, cc]);
        if size(bpoint, 1) >= 1 && bpoint(1) > 0 && isreal(bpoint(1)) && bpoint(1) < gamma_tr
            gamma(end+1, 1) = bpoint(1);
            IM(end+1, 1) = n;
        end
        if size(bpoint, 1) >= 2 && bpoint(2) > 0 && isreal(bpoint(2)) && bpoint(2) < gamma_tr
            gamma(end+1, 1) = bpoint(2);
            IM(end+1, 1) = n;
        end
    end
    [gamma, ind_temp] = sort(gamma);
    IM = IM(ind_temp);
    local_minima = [];
    local_minima_values = [];
    tl = 0;
    while ~isempty(IM)
        tu = gamma(1);
        n = IM(1);
        gamma = gamma(2:end);
        IM = IM(2:end);
        success = false;
        if g'*s0 > 0
            % Ascent
            0;
        elseif s0'*B*s0 < 0
            local_minima(end+1) = tu;
            success = true;
        else
            % Minimizing the quadratic model in direction d
            tau = -((s0'*B*s0)\(g'*s0));
            if tau < tl
                0;
            elseif tau < tu
                local_minima(end+1) = tau;
                success = true;
            else
                local_minima(end+1) = tu;
                success = true;
            end
        end
        if sign(tu*(s0'*cmodel(n).H*s0) + cmodel(n).g'*s0) > 0
            B = B + mu*cmodel(n).H;
            g = g + mu*cmodel(n).g;
            constraint_changed = n;
        else
            B = B - mu*cmodel(n).H;
            g = g - mu*cmodel(n).g;
            constraint_changed = -n;
        end
        tl = tu;
    end
    tu = gamma_tr;
    if g'*s0 > 0
        % Ascent
        0;
    elseif s0'*B*s0 < 0
        local_minima(end+1) = tu;
    else
        % Minimizing the quadratic model in direction d
        tau = -((s0'*B*s0)\(g'*s0));
        if tau < tl
            0;
        elseif tau < tu
            local_minima(end+1) = tau;
        else
            local_minima(end+1) = tu;
        end
    end
    n_minima = length(local_minima);
    pred = 0;
    for n = 1:n_minima
        s_n = local_minima(n)*s0;
        pred_n = predict_descent(fmodel, cmodel, s_n, mu, []);
        if pred_n > pred
            pred = pred_n;
            s = s_n;
        end
    end
    if pred > 0
        status = true;
    else
        % Rollback
        s = s0;
        pred = pred0;
        status = false;
    end
else
    s = s0;
    pred = pred0;
    status = true;
end
