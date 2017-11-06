function dp = l1_pseudo_gradient(gfx, mu, current_constraints, ind_eviolated)

dp = 0;
for n = ind_eviolated'
    dp = dp + mu*(current_constraints(n).g);
end
dp = dp + gfx;

end