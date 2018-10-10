function [c, g, H] = get_model_matrices(model, m, no_shift)
%GET_MODEL_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 3
        no_shift = false;
    end
    

    polynomial_shifted = model.modeling_polynomials{m + 1};


    if ~no_shift
        % Shift to TR center
        tr_center_shift = model.points_shifted(:, model.tr_center);
        polynomial = shift_polynomial(polynomial_shifted, tr_center_shift);
    end
    
    [c, g, H] = get_matrices(polynomial);

end

