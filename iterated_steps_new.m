function [s, pred] = iterated_steps_new(fmodel, cmodel, radius, mu, ...
                                        x0, Q, R, ind_qr, lb, ub)

tol1 = 0.95*radius;
tol2 = 0.5*radius;
tol3 = 0.05*tol2;
tol4 = 0.1*tol3;

dim = size(x0, 1);

pred = 0;
s = zeros(dim, 1);
for k = 1:dim

    pg = l1_pseudo_gradient_new(fmodel, cmodel, mu, s, ind_qr);

    h = null_space_conjugate_gradient(fmodel, cmodel, mu, Q, R, ...
                                      ind_qr, x0, -pg, radius, lb, ub, s);
    pred_h = predict_descent(fmodel, cmodel, h, mu, []);

    v = l1_range_step(fmodel, cmodel, Q, R, mu, h, ind_qr, radius, ...
                      x0, lb, ub);
    hv = (h + v);
    pred_hv = predict_descent(fmodel, cmodel, hv, mu, []);
    if pred_h < pred
        1;
    end
    if pred_hv >= pred_h
        s = hv;
        pred = pred_hv;
    else
        s = h;
        pred = pred_h;
    end
    if norm(s) - radius > sqrt(eps)
        1;
    end
    if norm(s) > tol1 ...
            || (norm(s) > tol2 && norm(hv) < tol3) ...
            || norm(hv) < tol4
        break
    end
end