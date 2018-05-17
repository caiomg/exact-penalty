function [point, value] = minimize_polynomial(polynomial, radius)


dimension = polynomial.dimension;
coefficients = polynomial.coefficients;

[c0, g, H] = coefficients_to_matrices(dimension, coefficients);

[point, value] = ms_step(H, g, radius);

if norm(point) > radius
    point = point*(radius/norm(point));
    value = c0 + g'*point + 0.5*(point'*H*point);
end

end