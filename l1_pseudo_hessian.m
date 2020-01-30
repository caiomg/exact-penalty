function Hp = l1_pseudo_hessian(fmodel, cmodel, mu, is_eactive, multipliers)


    n_constraints = length(cmodel);
    if numel(is_eactive) ~= n_constraints
        error('cmg:inconsistent_dimensions', 'Wrong number of elements');
    end
    if sum(is_eactive) ~= numel(multipliers)n
        error('cmg:inconsistent_dimensions', 'Wrong number of elements');
    end

    Hfx = fmodel.H;
    Hc = zeros(size(Hfx));
    Hm = zeros(size(Hfx));
    mult_i = 0;
    for n = 1:n_constraints
        if is_eactive(n)
            mult_i = mult_i + 1;
            Hm = Hm + multipliers(mult_i)*cmodel(n).H;
        elseif cmodel(n).c > 0
            Hc = Hc + cmodel(n).H;
        end
    end
    Hp = (Hfx + Hm) + mu*Hc;
end