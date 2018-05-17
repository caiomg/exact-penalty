function polynomial = choose_pivot_polynomial(polynomials)

% Just trying to choose a polynomial by some criterion
% For now, trying the polynomial with highest norm
n_polynomials = length(polynomials);
max_n = 0;
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
    if c_norm > max_n
        max_n = c_norm;
        chosen = np;
    end
end

polynomial = polynomials(chosen);


end