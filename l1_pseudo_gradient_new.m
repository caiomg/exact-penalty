function pg = l1_pseudo_gradient_new(fmodel, cmodel, mu, s, ind_eactive)

    if nargin < 4 || isempty(s)
        s = zeros(size(fmodel.g));
    end
    if nargin < 5
        ind_eactive = [];
    end
    
    n_constraints = length(cmodel);

    vg = 0;
    for n = 1:n_constraints
        c_value = cmodel(n).c + (cmodel(n).g + 0.5*(cmodel(n).H*s))'*s;
        if c_value > 0 && ~sum(n == ind_eactive)
            vg = vg + (cmodel(n).g + cmodel(n).H*s);
        end
    end
    gfx = fmodel.g + fmodel.H*s;
    pg = gfx + mu*vg;

end