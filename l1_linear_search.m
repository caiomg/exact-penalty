function x = l1_linear_search(pseudo_gradient, h, x0, ...
                              current_constraints, p, epsilon, delta)


% Question: do I use epsilon here?

a = pseudo_gradient'*h;
n_constraints = size(current_constraints, 1);
gamma = Inf(n_constraints, 1);
Ik = [];
for n = 1:n_constraints
    if abs(current_constraints(n).c) > epsilon
        gamma(n) = -current_constraints(n).c/(current_constraints(n).g'*h);
        if gamma(n) > 0
            Ik(end+1, 1) = n;
        end
    end
end

while true

    if isempty(Ik)
        % Linear approximation is unbounded
        % I will try one full Newton step
        min_gamma = 1;
        break
    end
    min_gamma = Inf;
    ind_gamma = 0;
    for n = Ik'
        if gamma(n) < min_gamma
            min_gamma = gamma(n);
            ind_gamma = n;
        end
    end

    g = current_constraints(ind_gamma).g;
    a = a + abs(g'*h);

    if a < 0
        Ik = Ik(Ik ~= ind_gamma);
    else
        break
    end

end

p0 = p(x0);
if p(x0 + min_gamma*h) < p0 - delta
  x = x0 + min_gamma*h;
else
    m = h'*pseudo_gradient;
    while p(x0 + min_gamma*h) > p0 + 0.5*m*min_gamma
        min_gamma = 0.5*min_gamma;
    end
    x = x0 + min_gamma*h;


end