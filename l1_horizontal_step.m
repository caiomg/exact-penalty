function h = l1_horizontal_step(fmodel, cmodel, mu, x0, ind_eactive, Q, R, radius, bl, bu, multipliers)

    if nargin < 11 || isempty(multipliers)
        multipliers = zeros(size(ind_eactive));
    end
    if isempty(bl)
        bl = -inf(size(x0));
    end
    if isempty(bu)
        bu = inf(size(x0));
    end
    if find(x0 > bu | x0 < bl, 1)
        error()
    end
    tol_radius = 1e-6;
    
    dimension = size(x0, 1);
    n_constraints = length(cmodel);
    r_columns = size(R, 2);
    N = Q(:, r_columns+1:end);    

    Bv = zeros(dimension);
    Bm = zeros(dimension);
    im = 1;
    for n = 1:n_constraints
        if sum(n == ind_eactive)
            Bm = Bm + multipliers(im)*cmodel(n).H;
            im = im + 1;
        elseif cmodel(n).c > 0
            Bv = Bv + mu*(cmodel(n).H);
        end
    end
    Ba = zeros(dimension);

    B = (Bv + fmodel.H) + Bm;
    g0 = l1_pseudo_gradient(fmodel.g, mu, cmodel, ind_eactive, true);
    g = g0;
    x = x0;

    s = x - x0;
    while true
       while true
           if isempty(N)
               break
           end
            g = g0 + B*s;
            % Evaluating matrices on reduced space
            Br = N'*B*N;
            gr = N'*g;
            % One shot of More-Sorensen algorithm
            s0 = ms_step(Br, gr, radius);
            d = N*s0;
            l_newly_active = (x == bl & d < 0);
            u_newly_active = (x == bu & d > 0);
            for k = 1:dimension
                if l_newly_active(k)
                    b = zeros(dimension, 1);
                    b(k) = bl(k);
                    [Q, R] = qrinsert(Q, R, 1, b);
                elseif u_newly_active(k)
                    b = zeros(dimension, 1);
                    b(k) = bu(k);
                    [Q, R] = qrinsert(Q, R, 1, b);
                end                    
            end
            if sum(l_newly_active | u_newly_active) > 0
                r_columns = size(R, 2);
                N = Q(:, r_columns+1:end);    
            else
                break
            end
        end
       if isempty(N) || norm(d) == 0
           break
       end        
        lower_breakpoints = (bl - x)./d;
        upper_breakpoints = (bu - x)./d;
        tr_breakpoint = roots([d'*d, 2*s'*d, s'*s - radius^2]);
        tr_breakpoint = min(tr_breakpoint(tr_breakpoint > 0));
        if isempty(tr_breakpoint) || ~isreal(tr_breakpoint)
            tr_breakpoint = inf;
        end
        
        l_breakpoints = (lower_breakpoints > 0);
        u_breakpoints = (upper_breakpoints > 0);
        bp = min([lower_breakpoints(l_breakpoints); upper_breakpoints(u_breakpoints); tr_breakpoint]);
        t = minimize_until_breakpoint(B, g, d, bp);

        x = x + t*d;
        l_newly_active = (t == lower_breakpoints & d < 0);
        u_newly_active = (t == upper_breakpoints & d > 0);
        x(l_newly_active) = bl(l_newly_active);
        x(u_newly_active) = bu(u_newly_active);

        s = x - x0;
        
        if t < bp || t == tr_breakpoint || norm(s) - radius >= tol_radius ...
                || norm(s) - norm(s0) >= tol_radius
            break
        end
    end
    s2 = remove_increasing_components(s, cmodel, ind_eactive);
    if norm(s - s2) > 0
        1;
        if norm(s - s2) > 1e-5
            1;
        end
    end
    
    [h2, pred, status] = line_search_full_domain(fmodel, cmodel, mu, s2, ...
                                                 radius);
    % I should also compare with cauchy step
    if ~status
        h = zeros(size(h2));
    else
        h = project_to_bounds(x0 + h2, bl, bu) - x0;
    end

end

