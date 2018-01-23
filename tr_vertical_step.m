function v = tr_vertical_step(p, x0, Q, R, c1, h, radius)

persistent success;
persistent fail;
if ~exist('success', 'var') || isempty(success)
    success = 0;
    fail = 0;
end
if isempty(c1)
    v = zeros(size(h));
else
    % Newton step
    v0 = -Q*(R'\c1);
    v = v0;
    if norm(h + v) > radius
        v = v*(radius/(norm(h+v) + norm(h)));
        tr_obj = @(x) quadratic(Q*(R*R')*Q', Q*(R*c1), 0, x);
        tr_constraint = @(x) quadratic(eye(length(h)), h, 0.5*(h'*h - radius^2), x);
        nlcon = @(x) constraints({tr_constraint}, {}, x, 1);
        fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                       'SpecifyObjectiveGradient', true);
        v1 = fmincon(tr_obj, v,[],[],[],[],[],[], nlcon, fmincon_options);
        v2 = solve_tr_problem(Q*(R*R')*Q', Q*(R*c1) - Q*(R*R')*Q'*h, radius) - h;
%         if tr_obj(v) < tr_obj(v2) - sqrt(eps)
%             1;
%         end
        v = v2;
    end
    [v, status] = simple_backtracking_line_search(p, h, v, 0, 0.70);
    if status == false
        v = zeros(size(h));
        fail = fail + 1;
    else
        success = success + 1;
    end
    
end
    
end