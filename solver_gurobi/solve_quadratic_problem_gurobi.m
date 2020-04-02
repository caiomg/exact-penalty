function [x, fval, status] = solve_quadratic_problem_gurobi(H, g, c, Aineq, bineq, Aeq, ...
                                            beq, lb, ub, x0)


    n_ineqs = size(bineq, 1);
    n_eqs = size(beq, 1);
    dim = size(g, 1);

    model.Q = 0.5*sparse(H);
    model.obj = g;
    model.objcon = c;
    
    if n_ineqs == 0 && n_eqs == 0
        A = sparse(zeros(0, dim));
    else
        A = sparse([Aineq;
                     Aeq]);
    end
    model.A = A;
    model.rhs = [bineq;
                   beq];
    model.sense = [repmat('<', n_ineqs, 1);
                   repmat('=', n_eqs, 1)];
    if isempty(lb)
        % We have to always set lb: Gurobi defaults to lb = 0
        lb = -inf(dim, 1);
    end
    model.lb = lb;
    if ~isempty(ub)
        model.ub = ub;
    end
    model.start = x0;
    
    params.OutputFlag = 0;
    params.NonConvex = 2;
    params.TimeLimit = 2*dim;
    
    result = gurobi(model, params);
    if isfield(result, 'x')
        x = result.x;
        fval = result.objval;
        if strcmp(result.status, 'OPTIMAL')
            status = 0;
        else
            status = -1;
        end
    else
        [x, fval, status] = solve_quadratic_problem_matlab(H, g, ...
                                                      c, Aineq, ...
                                                      bineq, Aeq, ...
                                                      beq, lb, ub, x0);
    end

end