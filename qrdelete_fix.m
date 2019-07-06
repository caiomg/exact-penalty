function [Q, R] = qrdelete_fix(Q, R, k, orient)
    
    Q0 = Q;
    R0 = R;
    if nargin < 4
        [Q, R] = qrdelete(Q, R, k);
    else
        [Q, R] = qrdelete(Q, R, k, orient);
    end
    if ~istriu(R) % FIXING MATLAB BUG
        mantain_cols = 1:size(R0, 2);
        mantain_cols(k) = [];
        A = Q0*R0(:, mantain_cols);
        [Q, R] = qr(A);
    end
    
end
