function pg = l1_pseudo_gradient_new(fmodel, cmodel, mu, s, ind_qr)

    dim = size(fmodel.g, 1);
    if nargin < 4 || isempty(s)
        s = zeros(dim, 1);
    end
    if nargin < 5
        ind_qr = [];
    end
    assert(numel(s) == dim);
    
    n_constraints = length(cmodel);

    vg = 0;
    for n = 1:n_constraints
        c_value = cmodel(n).c + (cmodel(n).g + 0.5*(cmodel(n).H*s))'*s;
        if c_value > 0 && ~sum(n == ind_qr)
            vg = vg + (cmodel(n).g + cmodel(n).H*s);
        end
    end
    gfx = fmodel.g + fmodel.H*s;
    pg = gfx + mu*vg;

end