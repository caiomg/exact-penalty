function pred = predict_descent_with_multipliers(fmodel, con_model, ...
                                                 s, mu, ind_qr, multipliers)

f_change = (0.5*s'*fmodel.H*s + fmodel.g'*s);
c_change = 0;
m_change = 0;
n_constraints = length(con_model);
for k = 1:n_constraints
    % Add change from ALL constraints
    this_change = (0.5*s'*con_model(k).H*s + con_model(k).g'*s);
    if con_model(k).c > 0
        if this_change + con_model(k).c < 0
            this_change = -con_model(k).c;
        end
        c_change = c_change + this_change;
    elseif con_model(k).c + this_change > 0
        this_change = con_model(k).c + this_change;
        c_change = c_change + this_change;
    end
    % From the active ones, add multiplier contribution over
    m_index = find(ind_qr == k, 1);
    if ~isempty(m_index)
        this_m_change = ...
            multipliers(m_index)*(0.5*s'*con_model(k).H*s + con_model(k).g'*s);
        m_change = m_change + this_m_change;
    end
end

% Transforming raw change into descent
pred = -(f_change + mu*c_change + m_change);

end