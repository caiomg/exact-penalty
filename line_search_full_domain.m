function [s, pred] = line_search_full_domain(fmodel, cmodel, mu, s0, radius)

if nargin < 5 || isempty(radius)
    radius = norm(s);
end

pred0 = predict_descent(fmodel, cmodel, s0, mu, []);
if pred0 < 0
    B = fmodel.H;
    g = fmodel.g;
    n_constraints = length(cmodel);
    active = [];
    for n = 1:n_constraints
        if cmodel(n).c > 0 || ...
               (cmodel(n).c ==0 && ...
                ((cmodel(n).g'*s0 + 0.5*(s0'*cmodel(n).H*s0) > 0) || ...
                (0.5*cmodel(n).g'*s0 + 0.125*(s0'*cmodel(n).H*s0) > 0)))
            g = g + mu*cmodel(n).g;
            B = B + mu*cmodel(n).H;
            active(end+1) = n;
        end
    end
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
    N = eye(max(size(s0)));
    for n = 1:n_constraints
        if cmodel(n).c ~= 0
            Hc = cmodel(n).H;
            gc = cmodel(n).g;
            cc = cmodel(n).c;
            bpoint = roots([0.5*(s0'*(N'*Hc*N)*s0),  gc'*N*s0, cc]);
            if size(bpoint, 1) >= 1 && bpoint(1) > 0 && isreal(bpoint(1)) && bpoint(1) < gamma_tr
                gamma(end+1, 1) = bpoint(1);
                IM(end+1, 1) = n;
            end
            if size(bpoint, 1) >= 2 && bpoint(2) > 0 && isreal(bpoint(2)) && bpoint(2) < gamma_tr
                gamma(end+1, 1) = bpoint(2);
                IM(end+1, 1) = n;
            end
        end
    end
    [gamma, ind_temp] = sort(gamma);
    IM = IM(ind_temp);

    while ~isempty(IM)
        tu = gamma(1);
        n = IM(1);
        gamma = gamma(2:end);
        IM = IM(2:end);
        w = zeros(size(s0));

        if g'*s0 > 0
            % Ascent
            break
        elseif s0'*B*s0 < 0
            segment_descent = (0.5*(s0'*B*s0)*alpha_s^2 + g'*s0*alpha_s) - (0.5*(s0'*B*s0)*tu^2 + g'*s0*tu);
            total_descent = total_descent + segment_descent;
            alpha_s = tu;
        else
            % Minimizing the quadratic model in direction d
            tau = -((s0'*B*s0)\(g'*s0));
            if tau < alpha_s
                % does this occur?
                s = w + alpha_s*s0;
                fs = total_descent;
                break;
            elseif tau < tu
                segment_descent = (0.5*(s0'*B*s0)*alpha_s^2 + g'*s0*alpha_s) - (0.5*(s0'*B*s0)*tau^2 + g'*s0*tau);
                total_descent = total_descent + segment_descent;
                alpha_s = tau;
                s = w + alpha_s*s0;
                fs = total_descent;
                % change the active set?
                break
            else
                segment_descent = (0.5*(s0'*B*s0)*alpha_s^2 + g'*s0*alpha_s) - (0.5*(s0'*B*s0)*tu^2 + g'*s0*tu);
                total_descent = total_descent + segment_descent;
                alpha_s = tu;
            end
        end
        if sign(alpha_s*((s0'*N')*cmodel(n).H*(N*s0)) + cmodel(n).g'*N*s0) > 0
            B = B + mu*N'*cmodel(n).H*N;
            g = g + mu*N'*cmodel(n).g;
            constraint_changed = n;
        else
            B = B - mu*N'*cmodel(n).H*N;
            g = g - mu*N'*cmodel(n).g;
            constraint_changed = -n;
        end
    end    
    if ~exist('s', 'var')
        if g'*s0 > 0
            % Ascent
            s = alpha_s*s0;
            fs = total_descent;
        elseif s0'*B*s0 < 0
            segment_descent = (0.5*(s0'*B*s0)*alpha_s^2 + g'*s0*alpha_s) - (0.5*(s0'*B*s0)*gamma_tr^2 + g'*s0*gamma_tr);
            alpha_s = gamma_tr;
            s = alpha_s*s0;
            total_descent = total_descent + segment_descent;
            fs = total_descent;
        else
            % Minimizing the quadratic model in direction d
            tau = -(g'*s0)/(s0'*B*s0);
            if tau < alpha_s
                % does this occur?
                s = alpha_s*s0;
                fs = total_descent;
            elseif tau < gamma_tr
                segment_descent = (0.5*(s0'*B*s0)*alpha_s^2 + g'*s0*alpha_s) - (0.5*(s0'*B*s0)*tau^2 + g'*s0*tau);
                total_descent = total_descent + segment_descent;
                alpha_s = tau;
                s = alpha_s*s0;
                fs = total_descent;
            else
                segment_descent = (0.5*(s0'*B*s0)*alpha_s^2 + g'*s0*alpha_s) - (0.5*(s0'*B*s0)*gamma_tr^2 + g'*s0*gamma_tr);
                total_descent = total_descent + segment_descent;
                alpha_s = gamma_tr;
                s = alpha_s*s0;
                fs = total_descent;
            end
        end
    end
    pred = fs;
else
    s = s0;
    pred = pred0;
end
