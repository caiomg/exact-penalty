function [s, pred] = iterated_steps_second_order(fmodel, cmodel, ...
                                                 radius, mu, x0, Q, ...
                                                 R, ind_qr, lb, ub)

    tol1 = 0.95*radius;
    tol2 = 0.5*radius;
    tol3 = 0.05*tol2;
    tol4 = 0.1*tol3;


    d = null_space_step_complete(fmodel, cmodel, mu, x0, ind_qr, Q, R, ...
                                 radius, lb, ub, zeros(size(ind_qr)));

    h = null_space_conjugate_gradient(fmodel, cmodel, mu, Q, R, ind_qr, ...
                                      x0, d, radius, lb, ub);
    pred_h = predict_descent(fmodel, cmodel, h, mu, []);

    v = l1_range_step(fmodel, cmodel, Q, R, mu, h, ind_qr, radius, x0, lb, ub);
    hv = (h + v);
    pred_hv = predict_descent(fmodel, cmodel, hv, mu, []);

    if pred_hv >= pred_h
        s = hv;
        pred = pred_hv;
    else
        s = h;
        pred = pred_h;
    end
end