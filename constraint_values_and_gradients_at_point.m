function [nphi, nA] = constraint_values_and_gradients_at_point(cmodel, ind, s)

    if isempty(ind)
        ind = 1:length(cmodel);
    end
    n_vars = size(s, 1);
    n_considered = sum(ind);
    nphi = zeros(n_considered, 1);
    nA = zeros(n_vars, n_considered);
    constraint_numbers = find(ind);
    for n = n_considered:-1:1
        nphi(n, 1) = cmodel(constraint_numbers(n)).c...
            + cmodel(constraint_numbers(n)).g'*s ...
            + 0.5*(s'*cmodel(constraint_numbers(n)).H*s);
        nA(:, n) = cmodel(constraint_numbers(n)).g ...
            + cmodel(constraint_numbers(n)).H*s;
    end

end
