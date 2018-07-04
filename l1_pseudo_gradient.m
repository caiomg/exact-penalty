function dp = l1_pseudo_gradient(gfx, mu, cmodel, ind_eactive, dummy)

if nargin < 5
    error();
end

n_constraints = length(cmodel);

dp = 0;
for n = 1:n_constraints
    if cmodel(n).c > 0 && ~sum(n == ind_eactive)
        dp = dp + mu*(cmodel(n).g);
    end
end
dp = dp + gfx;

end