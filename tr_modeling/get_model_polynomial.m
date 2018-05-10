function p = get_model_polynomial(model, m)
%GET_MODEL_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here


if m == 1
    poly = model.model_polynomial;
else
    poly = model.other_polynomials{m - 1};
end

[c, g, H] = coefficients_to_matrices(poly.dimension, poly.coefficients);

p.c = c;
p.g = g/model.scale_factor_x;
p.H = H/(model.scale_factor_x)^2;


end

