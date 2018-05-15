function shift_cmodel(cmodel, s)

    n_constraints = length(cmodel);
    [c, g] = update_constraint_information(cmodel, 1:n_constraints, s);
    for k = 1:n_constraints
        cmodel.c = c(k);
        cmodel.g = g(:, k);
    end
    
end