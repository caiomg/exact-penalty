function Hp = l1_hessian(fmodel, cmodel, mu, s)

    n_constraints = length(cmodel);

    Hfx = fmodel.H;
    Hc = zeros(size(Hfx));
    for n = 1:n_constraints
        if cmodel(n).c + cmodel(n).g'*s + 0.5*(s'*cmodel(n).H*s) > 0
            Hc = Hc + cmodel(n).H;
        end
    end
    Hp = Hfx + mu*Hc;

end