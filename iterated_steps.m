function s = iterated_steps(fmodel, cmodel, radius, mu, x0, Q, R, ind_qr, bl, bu)

tol1 = 0.95*radius;
tol2 = 0.5*radius;
tol3 = 0.05*tol2;
tol4 = 0.1*tol3;

dim = size(x0, 1);

s = zeros(dim, 1);
for k = 1:dim
    fmodel_d = shift_model(fmodel, s);
    for m = 1:length(cmodel)
        cmodel_d(m) = shift_model(cmodel(m), s);
    end
    x = project_to_bounds(x0 + s, bl, bu);

    h = null_space_step_cg(fmodel_d, cmodel_d, mu, x, ...
                            ind_qr, Q, R, radius - norm(s), bl, bu);
    pred_h = predict_descent(fmodel_d, cmodel_d, h, mu, []);

    v = l1_range_step(fmodel, cmodel, Q, R, mu, s + h, ind_qr, ...
                      radius, x0, bl, bu);
    sk = (h + v);
    pred = predict_descent(fmodel_d, cmodel_d, sk, mu, []);
    if pred < pred_h
        sk = h;
        pred = pred_h;
    end
    s = s + sk;
    if norm(s) - radius > sqrt(eps)
        1;
    end
    if norm(s) > tol1 ...
            || (norm(s) > tol2 && norm(sk) < tol3) ...
            || norm(sk) < tol4
        break
    end
end