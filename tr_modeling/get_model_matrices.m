function [c, g, H] = get_model_matrices(model, m)
%GET_MODEL_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here


if m == 0
    poly = model.model_polynomial;
else
    poly = model.other_polynomials{m};
end

[c, g, H] = coefficients_to_matrices(poly.dimension, poly.coefficients);


% Scale
g = g/model.scale_factor_x;
H = H/(model.scale_factor_x)^2;

end

