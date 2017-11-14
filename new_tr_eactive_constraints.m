function activated = new_tr_eactive_constraints(con_models, Ii, dx, epsilon)

n_constraints = length(Ii);
activated = zeros(0, 1);
for k = 1:n_constraints
    ck = con_models(Ii(k));
    value = 0.5*(dx'*ck.H*dx) + ck.g'*dx + ck.c;
    if abs(value) < epsilon
        activated(end+1, 1) = Ii(k);
    end
end
