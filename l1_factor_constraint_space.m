function [Q, R, ind_qr, lb_active, ub_active] = ...
    l1_factor_constraint_space(cmodel, con_eactive, lb_active, ub_active)
% L1_FACTOR_CONSTRAINT_SPACE - 
%   
    dim = size(cmodel(1).g, 1);
    A = [cmodel(con_eactive).g];
    I = eye(dim);
    A = [A, -I(:, lb_active), I(:, ub_active)];
    [Q, R] = qr(A);
    ind_qr = (1:length(cmodel))';
    ind_qr = ind_qr(con_eactive);
    
end

