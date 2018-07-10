function s = remove_increasing_components(s0, cmodel, ind_eactive)

n_cons = length(ind_eactive);
s = s0;
if n_cons ~= 0
    B = [cmodel(ind_eactive).g]';
    r_sum = zeros(n_cons, 1);
    r = max(0, B*s);
    while norm(r, inf) ~= 0
        lastwarn('');
        correction = - (B\r_sum);
        if ~isempty(lastwarn())
            1; % debuging
        end
        if norm(correction, inf) == 0
            break
        end
        s = s0 + correction;
        r = max(0, B*s);
        r_sum = r_sum + r;
        r_sum(r > 0) = 1.01*r_sum(r > 0);
    end
end


end