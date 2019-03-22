function Hp = l1_pseudo_hessian(fmodel, cmodel, mu, ind_qr, multipliers)


    n_constraints = length(cmodel);
    if length(ind_qr) ~= length(multipliers)
        1;
    end

    Hfx = fmodel.H;
    Hc = zeros(size(Hfx));
    Hm = zeros(size(Hfx));
    mult_i = 0;
    for n = 1:n_constraints
        if (isempty(find(ind_qr == n, 1)) ...
            && cmodel(n).c > 0)
            Hc = Hc + cmodel(n).H;
        elseif ~isempty(find(ind_qr == n, 1))
            mult_i = mult_i + 1;
            Hm = Hm + multipliers(mult_i)*cmodel(n).H;
        end
    end
    Hp = Hfx + mu*Hc + Hm;
end