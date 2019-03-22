function Hp = l1_hessian(fmodel, cmodel, mu, s, ind_eactive)

    if nargin < 5
        ind_eactive = [];
    end
    
    n_constraints = length(cmodel);

    Hfx = fmodel.H;
    Hc = zeros(size(Hfx));
    for n = 1:n_constraints
        
        if (isempty(find(ind_eactive == n, 1)) ...
            && cmodel(n).c + cmodel(n).g'*s + 0.5*(s'*cmodel(n).H*s)> 0)
            Hc = Hc + cmodel(n).H;
        end
    end
    Hp = Hfx + mu*Hc;

end