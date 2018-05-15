function v = tr_vertical_step_cg_origin(fmodel, cmodel, mu, h, ind_eactive, ind_eviolated, radius, Q, R)


opts_lt.LT = true;
tol_r = 1e-8;
dim = size(fmodel.g, 1);
n_cons = length(cmodel);
ls_steps = 10;
ls_factor = 0.75;
pv = @(s) -predict_descent(fmodel, cmodel, s, mu);

if ~isempty(ind_eactive)
    [nphi, Ah] = update_constraint_information(cmodel, 1:n_cons, h);
    Ah = Ah(:, ind_eactive);

    ind_pgrad = (1:n_cons)';
    ind_pgrad = ind_pgrad(nphi > 0);
    nphi = nphi(ind_eactive);
    [~, B] = update_constraint_information(cmodel, ind_pgrad, h);
    pp_grad = (fmodel.g + fmodel.H*h) + mu*sum(B, 2);


    gk0 = Ah*nphi;
    A = Q*R;
    gk = A*nphi;
    
%     if gk0'*pp_grad > 0 && gk0'*gk > 0
%         gk = gk0;
%     end
%     g_mn = A'\nphi;

    if gk'*pp_grad < 0
        % Try reseting the direction
        gk = gk - pp_grad*(gk'*pp_grad)/(pp_grad'*pp_grad);
%         gk = A'\nphi; 
%         if gk'*pp_grad < 0
%             gk = gk - pp_grad*(gk'*pp_grad)/(pp_grad'*pp_grad);
%         end
    end
    H = A*A';
    v = zeros(dim, 1);
    pk = -gk;
    % Tolerance!!
    for k = 1:dim
        if norm(h+v) - radius < tol_r
            break;
        end
        remaining_radius = roots([pk'*pk, 2*(h + v)'*pk, (h + v)'*(h + v) - radius^2]);
        radiusk = max(remaining_radius)*norm(pk);
        tpoint = pk*(radiusk/norm(pk));
        fmodelk = shift_model(fmodel, h + v);
        for n = 1:n_cons
            cmodelk(n) = shift_model(cmodel(n), h + v);
        end
        try
        [sk, ~, status] = line_search_full_domain(fmodelk, cmodelk, mu, tpoint, radiusk);
        catch erro
            rethrow(erro);
        end
        if status
            v = v + sk;
        end
        if ~status
           % Stop iterating
           break
        end
        gn = gk + (H*sk);
        b = (gn'*gn)/(gk'*gk);
        pk = -gn + b*pk;
        gk = gn;
    end

    if norm(pk, inf) == 0
        1;
    end

    if norm(pk, inf) ~= 0 && norm(v) > radius
        1;
%         % We have to truncate step at tr border
%         a = roots([pk'*pk, 2*v'*pk, v'*v - radius^2]);
%         a = max(a);
%         if isempty(a) || ~isreal(a) || a < 0
%            error();
%         end
%         v = v + a*pk;
    end
else
    v = zeros(size(h));
end




end



