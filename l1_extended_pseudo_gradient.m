function g = l1_extend_pseudo_gradient(mu, current_constraints, h, ind_eactive)

g = 0;
for n = ind_eactive'
    if current_constraints(n).g'*h < 0
        g = g - (current_constraints(n).g)/mu;
    end
end
g = g;

end