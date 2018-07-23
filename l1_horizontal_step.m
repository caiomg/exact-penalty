function [h, pred] = l1_horizontal_step(fmodel, cmodel, mu, x0, ind_eactive, Q, R, radius, bl, bu, multipliers)

    if nargin < 11 || isempty(multipliers)
        multipliers = zeros(size(ind_eactive));
        multipliers_provided = false;
    else
        multipliers_provided = true;
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
    tol_ort = 1e-5;
    
    dimension = size(x0, 1);
    n_constraints = length(cmodel);
    ind_csv = (1:n_constraints);
    ind_csv(ind_eactive) = [];

    r_columns = size(R, 2);
    N = Q(:, r_columns+1:end);    
    g0 = l1_pseudo_gradient(fmodel.g, mu, cmodel, ind_eactive, true);
    d0 = -(N*(N'*g0));

    Bv = zeros(dimension);
    Bm = zeros(dimension);
    im = 1;
    for n = 1:n_constraints
        if sum(n == ind_eactive)
            Bm = Bm + multipliers(im)*cmodel(n).H;
            im = im + 1;
        elseif cmodel(n).c > 0
            Bv = Bv + (cmodel(n).H);
        end
    end
    Ba = zeros(dimension);
    ba_included = false(size(ind_eactive, 1), 1);
%     if norm(multipliers, 'inf') == 0 && ~isempty(ind_eactive)
%         for n = ind_eactive'
%             if cmodel(n).c >= 0 && d0'*cmodel(n).H*d0 > 0
%                 Ba = Ba + (cmodel(n).H);
%             end
%         end
%     end

    while true
        B = (fmodel.H + mu*Bv) + Bm + mu*Ba;
        g = g0;
        x = x0;

        s = zeros(size(x0));
        l_bounds_included = false(dimension, 1);
        u_bounds_included = false(dimension, 1);
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
               inserted = false;
               for k = 1:dimension
                   if l_newly_active(k)
                       b = zeros(dimension, 1);
                       b(k) = 1;
                       norm_b = norm((Q*R)*(R\(Q'*b)) - b, 1);
                       if norm_b > tol_ort
                           [Q, R] = qrinsert(Q, R, 1, b);
                           inserted = true;
                           l_bounds_included(k) = true;
                           break
                       end
                   elseif u_newly_active(k)
                       b = zeros(dimension, 1);
                       b(k) = -1;
                       norm_b = norm((Q*R)*(R\(Q'*b)) - b, 1);
                       if norm_b > tol_ort
                           [Q, R] = qrinsert(Q, R, 1, b);
                           inserted = true;
                           u_bounds_included(k) = true;
                           break
                       end
                   end                    
                end
                if inserted
                    r_columns = size(R, 2);
                    N = Q(:, r_columns+1:end);
                    N(l_bounds_included, :) = zeros(sum(l_bounds_included), size(N, 2));
                    N(u_bounds_included, :) = zeros(sum(u_bounds_included), size(N, 2));
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

            l_newly_active = (t == lower_breakpoints & d < 0);
            u_newly_active = (t == upper_breakpoints & d > 0);
            s = correct_step_to_bounds(x0, s + t*d, bl, bu, l_newly_active, u_newly_active);

            x = x0 + s;
            x(l_bounds_included) = bl(l_bounds_included);
            x(l_newly_active) = bl(l_newly_active);
            x(u_bounds_included) = bu(u_bounds_included);
            x(u_newly_active) = bu(u_newly_active);

            if t < bp || t == tr_breakpoint || norm(s) - radius >= tol_radius ...
                    || norm(s) - norm(s0) >= tol_radius
                break
            end
        end

%         s2 = remove_increasing_components(s, cmodel, ind_eactive);


        pred_h = predict_descent(fmodel, cmodel, s, mu, ind_eactive);
        if ~multipliers_provided || radius - norm(s) < tol_radius
            [h, pred, status] = line_search_full_domain(fmodel, cmodel, mu, s, ...
                                                         radius);
            if ~status
                h = zeros(size(h));
            end
            break
        else
            h = s;
            pred = predict_descent(fmodel, cmodel, h, mu, []);
            break
        end
    end

    if ~find(x < bl | x > bu, 1)
        % If x inside bounds we can try adjusting the step
        h = correct_step_to_bounds(x0, h, bl, bu);
    end
    



end

