function dp = l1_pseudo_gradient(gfx, mu, current_constraints, ind_eviolated)

dp = 0;
for n = ind_eviolated'
    dp = dp - (current_constraints(n).g)/mu;
end
dp = dp + gfx;

end