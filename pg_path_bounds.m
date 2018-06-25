function [s, fs] = pg_path_bounds(H, g, x0, s0, bl, bu, radius)

if isempty(bl)
    bl = -inf(size(x));
end
if isempty(bu)
    bu = inf(size(x));
end

dim = size(g, 1);
tol = dim*10*eps;

x = x0 + s0;

lower_hits = x < bl;
upper_hits = x > bu;
while sum(lower_hits | upper_hits)
    not_hits = ~(lower_hits | upper_hits);
    x = min(max(x, bl), bu);
    gs = g + H*x;
    g_red = gs(not_hits);
    try
        H_red = H(not_hits, not_hits);
    catch err
        rethrow(err)
    end
    d = s0(not_hits);
    if isempty(d) || norm(d) < tol || norm(x) >= radius || norm(g_red) < tol
        break
    end
    
    lower_breakpoints = (bl(not_hits) - x(not_hits))./d;
    next_l_ind= (1:dim)';
    next_l_ind = next_l_ind(not_hits);
    next_l_ind = next_l_ind(lower_breakpoints > 0);
    lower_breakpoints = lower_breakpoints(lower_breakpoints > 0);
    [next_l_bp, ind] = min(lower_breakpoints);
    next_l_ind = next_l_ind(ind);
    
    upper_breakpoints = (bu(not_hits) - x(not_hits))./d;
    next_u_ind= (1:dim)';
    next_u_ind = next_u_ind(not_hits);
    next_u_ind = next_u_ind(upper_breakpoints > 0);
    upper_breakpoints = upper_breakpoints(upper_breakpoints > 0);
    [next_u_bp, ind] = min(upper_breakpoints);
    next_u_ind = next_u_ind(ind);

    tr_breakpoint = roots([d'*d, 2*x(not_hits)'*d, x'*x - radius^2]);
    if ~isreal(tr_breakpoint)
        tr_breakpoint = inf;
    end
    tr_breakpoint = min(tr_breakpoint(tr_breakpoint > 0));

    t = -g_red'*d/(d'*H_red*d);
    t = min([t, next_l_bp, next_u_bp, tr_breakpoint]);
    if t < 0
        break
    end
    x(not_hits) = x(not_hits) + t*d;
    if t == tr_breakpoint
        break
    elseif (isempty(next_l_bp) || t ~= next_l_bp) && (isempty(next_u_bp) || t ~= next_u_bp)
        break
    else
        if t == next_l_bp
            lower_hits(next_l_ind) = true;
            x(next_l_ind) = bl(next_l_ind); % small correction
        end
        if t == next_u_bp
            upper_hits(next_u_ind) = true;
            x(next_u_ind) = bu(next_u_ind); % small correction
        end
    end
end

s = x - x0;
fs = g'*s + 0.5*(s'*H*s);
if sum(x > bu | x < bl)
    error('cmg:bounds_violated', 'Bounds not satisfied');
end

end