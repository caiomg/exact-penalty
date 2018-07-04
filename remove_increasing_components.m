function s = remove_increasing_components(s0, cmodel, ind_eactive)

n_cons = length(ind_eactive);
if n_cons == 0
    s = s0;
else
    B = [cmodel(ind_eactive).g]';
    r_sum = zeros(n_cons, 1);
    while true
        s = s0 - (B\r_sum);
        r = max(0, B*s);
        r_sum = r_sum + r;
        r_sum(r > 0) = 1.01*r_sum(r > 0);
        if norm(r, inf) == 0
            break
        end
    end
end


end