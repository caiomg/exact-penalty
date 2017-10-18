function x = l1_linear_search(B, pseudo_gradient, h, x0, ...
                              current_constraints, p, mu, epsilon, ...
                              delta, ind_eactive)


% Question: do I use epsilon here?

a = pseudo_gradient'*h;
a2 = h'*B*h;

n_constraints = size(current_constraints, 1);
gamma = Inf(n_constraints, 1);
Ik = [];
for n = 1:n_constraints
    if ~isempty(find(ind_eactive == n, 1))
        % Not testing constraints which are already active
        continue
    end
    if current_constraints(n).c ~= 0
        Hc = current_constraints(n).H;
        gc = current_constraints(n).g;
        cc = current_constraints(n).c;
        bpoint = roots([0.5*(h'*Hc*h), gc'*h, cc]);
        % Minimum of the roots
        gamma_n = min(bpoint(bpoint > 0));
        if ~isempty(gamma_n) && isreal(gamma_n)
            gamma(n) = gamma_n;
            Ik(end+1, 1) = n;
        end
    end
end
if isempty(Ik)
    alpha = 1;
else
    alpha = 0;
end

while ~isempty(Ik)

    [min_gamma, ind_temp] = min(gamma(Ik));
    ind_gamma = Ik(ind_temp);
    Ik = Ik(Ik ~= ind_gamma);
    if a > 0
        break
    elseif a2 <= 0
        % Go until the end
        alpha = min_gamma;
    else
        % Find minimum
        tau = -a/a2;
        if tau < alpha
            break
        elseif tau < min_gamma
            alpha = tau;
            break
        else
            alpha = min_gamma;
        end
    end
    g = current_constraints(ind_gamma).g;
    sigma = sign(g'*h);
    a = a + sigma*(g'*h)/mu;
    Hc = current_constraints(ind_gamma).H;
    a2 = a2 + sigma*(h'*Hc*h)/mu;
end

x = x0 + alpha*h;

% g = pseudo_gradient + l1_extended_pseudo_gradient(mu, current_constraints, ...
%                                                   h, ind_eactive);
% m = h'*g;
% p0 = p(x0);
% if p(x0 + alpha*h) < p0 - delta
%   x = x0 + alpha*h;
% else
%     while p(x0 + alpha*h) > p0 - 0.5*m^2
%         if norm((x0 + (0.5*alpha)*h) - x0, 'inf') == 0
%             break
%         end
%         alpha = 0.5*alpha;
%     end
%     x = x0 + alpha*h;
% end
