function v = tr_vertical_step(Q, R, c1, h, radius)

% Newton step
v = -Q*(R'\c1);

if norm(h + v) > radius
    v = v*(radius/(norm(h+v) + norm(h)));
    tr_obj = @(x) quadratic(Q*(R*R')*Q', Q*(R*c1), 0, x);
    tr_constraint = @(x) quadratic(Q*(R*R')*Q', Q*(R*c1), 0.5*(h'*h - radius^2), x);
    nlcon = @(x) constraints({tr_constraint}, {}, x, 1);
    fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                   'SpecifyObjectiveGradient', true);
    v = fmincon(tr_obj, v,[],[],[],[],[],[], nlcon, fmincon_options);
end
    
    
end