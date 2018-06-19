function [point, value] = minimize_polynomial(polynomial, radius)


dimension = polynomial.dimension;
coefficients = polynomial.coefficients;

[c0, g, H] = coefficients_to_matrices(dimension, coefficients);

[point, value] = ms_step(H, g, radius);
npoint = norm(point);

if npoint > radius
    point = (point/npoint)*radius;
    value = c0 + g'*point + 0.5*(point'*H*point);
end

end