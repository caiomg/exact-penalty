function [degenerate, d, con_confirmed, lb_confirmed, ub_confirmed] = ...
        l1_test_degeneracy(p_grad, Q_ext, R_ext, con_eactive, lb_active, ub_active, mu)
% L1_TEST_DEGENERACY - 
%   
    tol_independence = sqrt(eps(1));
    tol_activation = 1e-9;

    [dim, r_cols] = size(R_ext);
    if r_cols == 0 ...
            || (r_cols == 1 && norm(R_ext)> tol_independence) ...
            || (r_cols <= dim && min(abs(diag(R_ext))) > tol_independence)
        % Not degenerate
        degenerate = false;
        d = [];
        con_confirmed = con_eactive;
        lb_confirmed = lb_active;
        ub_confirmed = ub_active;
    else
        degenerate = true;
        lb = zeros(dim, 1);
        lb(~lb_active) = -ones(sum(~lb_active), 1);
        ub = zeros(dim, 1);
        ub(~ub_active) = ones(sum(~ub_active), 1);
        n_active_bounds = sum(lb_active | ub_active);
        n_constraints = r_cols - n_active_bounds;
        A_con = mu*(Q_ext*R_ext(:,1:n_constraints))';
        
         [d, lambda] = solve_linear_problem(p_grad, A_con, zeros(n_constraints, 1), [], [], lb, ub);
         %d = linprog(p_grad, A_con, zeros(n_constraints, 1), [], [], lb, ub);
        
%         intlinprog_problem.f = [p_grad; 0];
%         intlinprog_problem.solver = 'intlinprog';
%         intlinprog_problem.Aineq = [A_con, zeros(size(A_con, 1), 1)];
%         intlinprog_problem.bineq = zeros(n_constraints, 1);
%         intlinprog_problem.lb = [lb; 0];
%         intlinprog_problem.ub = [ub; 0];
%         intlinprog_problem.intcon = size(p_grad) + 1;
%         intlinprog_problem.options = optimoptions('intlinprog', ...
%                                                   'Display', 'off');
%         %, 'Algorithm', 'dual-simplex');
%         [d, ~, exitflag, ~] = intlinprog(intlinprog_problem);
%         % con_confirmed(con_eactive) = lambda.ineqlin > 0;
%         % lb_confirmed = lb_active & lambda.lower(1:end-1) > 0;
%         % ub_confirmed = ub_active & lambda.upper(1:end-1) > 0;
%         d = d(1:end-1);
        con_confirmed = con_eactive;
        con_confirmed(con_eactive) = A_con*d > 0;
        lb_confirmed = lb_active & (d < -1e-6);
        ub_confirmed = ub_active & (d > 1e-6);
    end
        
end