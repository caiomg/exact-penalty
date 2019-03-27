function constraints_models = ...
        extract_constraints_from_tr_model(model, con_lb, con_ub)


n_constraints = size(model.fvalues, 1) - 1;
assert(n_constraints == size(con_lb, 1));
assert(n_constraints == size(con_ub, 1));


constraints_models = [];
ind = 0;
for m = 1:n_constraints
    if isfinite(con_ub(m))
        ind = ind+1;
        [constraints_models(ind, 1).c, ...
         constraints_models(ind, 1).g, ...
         constraints_models(ind, 1).H] = get_model_matrices(model, m);
        constraints_models(ind, 1).c = ...
            constraints_models(ind).c - con_ub(m);
    end
    if isfinite(con_lb(m))
        ind = ind+1;
        [constraints_models(ind, 1).c, ...
         constraints_models(ind, 1).g, ...
         constraints_models(ind, 1).H] = get_model_matrices(model, m);
        constraints_models(ind).c = con_lb(m) - constraints_models(ind).c;
        constraints_models(ind).g = -constraints_models(ind).g;
        constraints_models(ind).H = -constraints_models(ind).H;
    end
end


end

