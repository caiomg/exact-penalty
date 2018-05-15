function v = tr_vertical_step(fmodel, cmodel, mu, h, ind_eactive, ind_eviolated, radius, Q, R)


opts_lt.LT = true;
tol_r = 1e-8;
dim = size(fmodel.g, 1);
n_cons = length(cmodel);
ls_steps = 10;
ls_factor = 0.75;
pv = @(s) -predict_descent(fmodel, cmodel, s, mu);

if ~isempty(ind_eactive)
    % vertical step without TR (Newton step)
    nphi = [cmodel(ind_eactive).c]';
    vh = linsolve(R', -nphi, opts_lt);
    v1 = Q*vh;
    % steepest descent direction for v step
    v2 = Q*(R*nphi);
    
    r1 = max(roots([v1'*v1, 2*h'*v1, h'*h - 1.5*radius^2]));
    r2 = max(roots([v2'*v2, 2*h'*v2, h'*h - 1.5*radius^2]));
    if isempty(r1) || ~isreal(r1)
        v1 = [];
        r1 = 0;
    end
    if isempty(r2) || ~isreal(r2)
        v2 = [];
        r2 = 0;
    end
    S = diag([r1, r2]);

    V = [v1, v2]*S;
    l_side = V'*Q*R;
    
    try
    g_bid = l_side*nphi;
    H_bid = l_side*l_side';
    catch erro
        1;
    end

    % More-Sorensen step
    try
    va = solve_tr_problem(H_bid, g_bid, 1);
    catch erro
        1;
    end
    v_before = V*va;
    
    fmodelk = shift_model(fmodel, h);
    for n = 1:n_cons
        cmodelk(n) = shift_model(cmodel(n), h);
    end
    r = max(roots([v_before'*v_before, 2*h'*v_before, h'*h - radius^2]));
    if isempty(r) || ~isreal(r)
        1;
    end
    radius_ls = r*norm(v_before);
    [v_after_ls, predv, status] = line_search_full_domain(fmodelk, cmodelk, mu, v_before, radius_ls);
    if status
        v = v_after_ls;
    else
        v = zeros(dim, 1);
    end
else
    v = zeros(size(h));
end




end



