function h = tr_horizontal_cg_step(fmodel, cmodel, mu, ind_eactive, ind_eviolated, radius, N, multipliers, ind_qr)


opts_lt.LT = true;
tol_r = 1e-8;
dim = size(fmodel.g, 1);
n_cons = length(cmodel);

Bv = zeros(dim);
for n = ind_eviolated'
    Bv = Bv + mu*(cmodel(n).H);
end
Bm = zeros(dim);
for n = 1:length(multipliers)
    Bm = Bm + multipliers(n)*(cmodel(ind_qr(n)).H);
end
H = Bv + Bm + fmodel.H;
pg = l1_pseudo_gradient(fmodel.g, mu, cmodel, ind_eviolated);

H_red = N'*H*N;
g_red = N'*pg;

u = zeros(size(g_red));
h = zeros(dim, 1);
gk = g_red;
pk = -gk;
for k = 1:dim
    if norm(u) >= radius - tol_r
        break
    end
    remaining_radius = roots([pk'*pk, 2*u'*pk, u'*u - radius^2]);
    radiusk = max(remaining_radius)*norm(pk);
    tpoint = (pk*(radiusk/norm(pk)));
    fmodelk = shift_model(fmodel, h);
    for n = 1:n_cons
        cmodelk(n) = shift_model(cmodel(n), h);
    end

    try
    [sk, ~, status] = line_search_full_domain(fmodelk, cmodelk, mu, N*tpoint, radiusk);
    catch erro
        1;
    end
    if status
        u = u + N'*sk;
        h = h + sk;
    end
    if ~status
       % Stop iterating
       break
    end
    gn = gk + (H_red*(N'*sk));
    b = (gn'*gn)/(gk'*gk);
    pk = -gn + b*pk;
    gk = gn;
end

if norm(pk, inf) == 0
    1;
end

end



