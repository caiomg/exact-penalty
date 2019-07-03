function [Q, R] = qrdelete_fix(Q, R, k, orient)
    
    if nargin < 4
        [Q, R] = qrdelete(Q, R, k);
    else
        [Q, R] = qrdelete(Q, R, k, orient);
    end
    if ~istriu(R) % FIXING MATLAB BUG
       [Q, R] = qr(Q*R);
    end
    
end
