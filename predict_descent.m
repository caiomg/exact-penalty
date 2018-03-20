function pred = predict_descent(fmodel, con_model, s, mu, ind_eactive)

if nargin < 5
    ind_eactive = [];
end

f_change = (0.5*(s'*fmodel.H*s) + fmodel.g'*s);
c_change = 0;
n_constraints = length(con_model);
for k = 1:n_constraints
    if isempty(find(ind_eactive == k, 1))
        this_change = (0.5*(s'*con_model(k).H*s) + con_model(k).g'*s);
        if con_model(k).c > 0
            if this_change + con_model(k).c < 0
                this_change = -con_model(k).c;
            end
            c_change = c_change + this_change;
        elseif con_model(k).c + this_change > 0
            this_change = con_model(k).c + this_change;
            c_change = c_change + this_change;
        end
    end
end

pred = -(f_change + mu*c_change);

end