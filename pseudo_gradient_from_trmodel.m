function pg = pseudo_gradient_from_trmodel(model, mu)

    [~, g] = get_model_matrices(model, 0);
    current_constraints = extract_constraints_from_tr_model(model);
    gc_sum = zeros(size(g));
    n_constraints = length(current_constraints);
    for k = 1:n_constraints
        if current_constraints(k).c > 0
            gc_sum = gc_sum + current_constraints(k).g;
        end
    end
    pg = g + mu*gc_sum;

end