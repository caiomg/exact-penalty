function [nphi, nA] = update_constraint_information(cmodel, ind, s)

    n_vars = size(s, 1);
    n_considered = length(ind);
    nphi = zeros(n_considered, 1);
    nA = zeros(n_vars, n_considered);
    for n = n_considered:-1:1
        nphi(n, 1) = cmodel(ind(n)).c + cmodel(ind(n)).g'*s + 0.5*(s'*cmodel(ind(n)).H*s);
        nA(:, n) = cmodel(ind(n)).g + cmodel(ind(n)).H*s;
    end

end