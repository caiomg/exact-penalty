function v = tr_vertical_step_new(fmodel, cmodel, mu, h, ind_eactive, ind_eviolated, radius)



n_cons = length(cmodel);
pv = @(s) -predict_descent(fmodel, cmodel, s, mu);

if ~isempty(ind_eactive)
    [nphi, A] = update_constraint_information(cmodel, ind_eactive, h);
    [phiB, B] = update_constraint_information(cmodel, [ind_eactive(nphi > 0); ind_eviolated], h);
    pp_grad = fmodel.g + mu*sum(B, 2);


    vc = -(A*nphi);
    scale_vc = (vc'*vc)/((vc'*A)*(A'*vc));
    vc = scale_vc*vc;
    if norm(h + vc) > radius
        if radius - norm(h) > 100*eps
           scale_vc = roots([vc'*vc, 2*h'*vc, h'*h - radius^2]);
           scale_vc = max(scale_vc);
           if scale_vc < 0 || ~isreal(scale_vc)
               error();
           end
        else
            scale_vc = 0;
        end
       vc = vc*scale_vc;
    end

    Z = null(A');
    vb1 = A'\-nphi;
    
    if vb1'*pp_grad <= 0
        vb2 = zeros(size(vb1));
    else
        vb2 = Z*((pp_grad'*Z)\(-1.1*vb1'*pp_grad));
    end
    vb = vb1 + vb2;

    if norm(h + vb) <= radius
        [vi, status] = simple_backtracking_line_search(pv, h, vb, 0, 0.70);
        if status
            v = vi;
        else
            if vc'*pp_grad <= 0
                [vi, status] = simple_backtracking_line_search(pv, h, vc, 0, 0.70);
                if status
                    v = vi;
                else
                    v = zeros(size(vi));
                end
            else
                v = zeros(size(h));
            end
        end
    else
        if vc'*pp_grad <= 0
            [vc, status] = simple_backtracking_line_search(pv, h, vc, 0, 0.70);
            if ~status
                vc = zeros(size(vc));
            end
        else
            vc = zeros(size(vc));
        end
        h2 = h + vc;
        d2 = vb - h2;
        if radius - norm(h2) > 100*eps
            scale_d2 = roots([d2'*d2, 2*h2'*d2, h2'*h2 - radius^2]);
            scale_d2 = max(scale_d2);
            if scale_d2 < 0 || ~isreal(scale_d2)
                error();
            end
        else
            scale_d2 = 0;
        end
        d2 = d2*scale_d2;
        [di, status] = simple_backtracking_line_search(pv, h2, d2, 0, 0.70);
        if status
            v = vc + di;
        else
            v = vc;
        end
    end
        
%     if norm(h + vc) > radius
%         h1 = h + vc;
%         h1 = h1*(radius/norm(h1));
%         vc = h1 - h;
%         [vc, status] = simple_backtracking_line_search(pv, h, vc, 0, 0.70);
%         if ~status
%             vc = zeros(size(vc));
%         end
%         v = vc;
%     else
%         % Computing pseudo-gradient on current point
%         [~, B] = update_constraint_information(cmodel, ind_eviolated, h+vc);
%         pp_grad = fmodel.g + mu*sum(B, 2);
% 
%         Z = null(A');
%         w = Z\-pp_grad;
%         vb1 = A'\-nphi;
%         vb2 = Z*w;
% 
%         h2 = h + vb1;
%         if norm(h2) >= radius
%             h2 = h2*(radius/norm(h2));
%             vb = h2 - h;
%         else
%             scale_n = roots([vb2'*vb2, 2*vb2'*h2, h2'*h2 - radius^2]);
%             scale_n = max(scale_n);
%             if scale_n < 0 || ~isreal(scale_n)
%                 error();
%             end
%             vb = vb1 + scale_n*vb2;
%         end
% 
%         [vi, status] = simple_backtracking_line_search(pv, h + vc, vb - vc, 0, 0.70);
%         if status
%             v = vc + vi;
%         else
%             v = vc;
%         end
%     end
%     v2 = solve_tr_problem(A*A', A*nphi, radius);
%     if norm(h + v2) > radius
%         hv2 = (h + v2)*(radius/norm(h + v2));
%         v2 = hv2 - h;
%     end
%     [v, status] = simple_backtracking_line_search(pv, h, v2, 0, 0.70);
%     if ~status
%         v = zeros(size(v));
%     end

else
    v = zeros(size(h));
end




end



