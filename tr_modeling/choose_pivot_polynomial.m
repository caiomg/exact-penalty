function [polynomial, index] = choose_pivot_polynomial(polynomials)

% Just trying to choose a polynomial by some criterion
n_polynomials = length(polynomials);
min_n = inf;
chosen = 0;
for np = 1:n_polynomials
    dimension = polynomials(np).dimension;
    if norm(polynomials(np).coefficients(2:dimension+1), inf) >= 1
        % Preference for 'linear' polynomials
        chosen = np;
        break
    end
    c_norm = norm(polynomials(np).coefficients, 1);
    % Trying to find the polynomial with highest norm
    if c_norm < min_n
        min_n = c_norm;
        chosen = np;
    end
end

polynomial = polynomials(chosen);
index = chosen;

end