function pg = l1_pseudo_gradient(gfx, mu, cmodel, ind_eactive, dummy)

if nargin < 5
    error();
end

n_constraints = length(cmodel);

vg = 0;
for n = 1:n_constraints
    if cmodel(n).c > 0 && ~sum(n == ind_eactive)
        vg = vg + cmodel(n).g;
    end
end
pg = gfx + mu*vg;

end