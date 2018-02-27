function constraints_models = extract_constraints_from_tr_model(model)
%EXTRACT_CONSTRAINTS_FROM_TR_MODEL Summary of this function goes here
%   Detailed explanation goes here


n_constraints = size(model.fvalues, 1) - 1;
constraints_models = [];
for m = 1:n_constraints
    [constraints_models(m, 1).c, constraints_models(m, 1).g, ...
     constraints_models(m, 1).H] = get_model_matrices(model, m);
end


end

