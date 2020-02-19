function [nphi, nA] = constraint_values_and_gradients_at_point(cmodel, s)

    n_vars = size(s, 1);
    n_considered = numel(cmodel);
    nphi = zeros(n_considered, 1);
    nA = zeros(n_vars, n_considered);
    for n = 1:n_considered
        nphi(n, 1) = cmodel(n).c + cmodel(n).g'*s + 0.5*(s'*cmodel(n).H*s);
        nA(:, n) = cmodel(n).g + cmodel(n).H*s;
    end

end
