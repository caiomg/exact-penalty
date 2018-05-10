function [point, value] = minimize_polynomial(polynomial, radius)


dimension = polynomial.dimension;
coefficients = polynomial.coefficients;

[c0, g, H] = coefficients_to_matrices(dimension, coefficients);

[point, value] = ms_step(H, g, radius);

end