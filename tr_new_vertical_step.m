function v = tr_new_vertical_step(p, cmodel, h, radius, ind_constraints, p_gradient)

max_iter = 10; % backtracking iterations
min_dec = 0; % minimum decrease accepted
n_constraints = length(cmodel);

if ~isempty(ind_constraints)
    g = zeros(size(cmodel(1).g));
    B = zeros(size(cmodel(1).H));
    

    for m = ind_constraints'
        % New values for the model
        c1 = cmodel(m).c + cmodel(m).g'*h + 0.5*(h'*cmodel(m).H*h);
        g1 = (cmodel(m).g + cmodel(m).H*h);
        % Summing
        g = g + g1*c1;
        B = B + ((g1*g1') + cmodel(m).H*c1);
    end
    
    v = zeros(size(h));
    for n = 1:max_iter
        s = solve_tr_problem(B, g, radius);
        if p(s) <= p(h) - min_dec
            v = s - h;
            break
        else
            radius = radius/2;
        end
    end
else
    v = zeros(size(cmodel(1).g));    
end
end