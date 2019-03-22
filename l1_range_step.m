function v = l1_range_step(fmodel, cmodel, Q, R, mu, h, ind_qr, ...
                           radius, x0, lb, ub)
    
    % Actual calculation of the whole step
    v = l1_range_step_complete(fmodel, cmodel, Q, R, mu, h, ind_qr, ...
                               radius, x0, lb, ub);

    % Backtracking
    pred_h = predict_descent(fmodel, cmodel, h, mu);
    h_fmodel = shift_model(fmodel, h);
    for m = 1:length(cmodel)
        h_cmodel(m) = shift_model(cmodel(m), h);
    end

    tol_v = 100*eps(1);
    pred = predict_descent(h_fmodel, h_cmodel, v, mu);
    while norm(v, inf) > tol_v && pred < -0.5*pred_h
        v = 0.5*v;
        pred = predict_descent(h_fmodel, h_cmodel, v, mu);
    end
end
