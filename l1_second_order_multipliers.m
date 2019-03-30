function eta = l1_second_order_multipliers(fmodel, cmodel, pgradient, ...
                                             mu, Q, R, ind_qr, multipliers)
    
    r_cols = size(R, 2);
    R1 = R(1:r_cols, 1:r_cols);
    % 'Vertical step'
    v = Q*(R1\(-[cmodel(ind_qr).c]));
    
    % Second-order estimates
    W = l1_pseudo_hessian(fmodel, cmodel, mu, ind_qr, multipliers);
    Q1 = Q(1:r_cols, :);
    eta = -(R1\(Q1*(pgradient + W*v)));
    
    
end
    
    