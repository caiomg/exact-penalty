function [x, lambda] = solve_linear_problem(f, Aineq, bineq, Aeq, beq, lb, ub, x0)
    
    
    dim = size(f, 1);
    if nargin < 8 || isempty(x0)
        x0 = zeros(dim, 1);
    end
    if nargin < 7
        ub = [];
    end
    if nargin < 6
        lb = [];
    end
    if nargin < 4
        Aeq = [];
        beq = [];
    end
    if nargin < 2
        Aineq = [];
        bineq = [];
    end
    
    f_ipopt.objective = @(x) f'*x;
    f_ipopt.gradient = @(x) f;
    f_ipopt.hessian = @(x, sigma, lambda) sparse(zeros(dim));
    f_ipopt.hessianstructure = @() sparse(zeros(dim));

    ipopt_options.ipopt.hessian_constant = 'yes';
    if ~isempty(lb)
        ipopt_options.lb = lb;
    end
    if ~isempty(ub)
        ipopt_options.ub = ub;
    end
    if ~isempty(Aineq)
        n_ineq = numel(bineq);
        ipopt_options.cu = bineq;
        ipopt_options.cl = -inf(n_ineq, 1);
        f_ipopt.constraints = @(x) Aineq*x;
        f_ipopt.jacobian = @(x) sparse(Aineq);
        f_ipopt.jacobianstructure = @(x) sparse(ones(n_ineq, dim));
        ipopt_options.ipopt.jac_d_constant = 'yes';
    end
    if ~isempty(Aeq)
        ipopt_options.ipopt.jac_c_constant = 'yes';
        error('Not implemented');
    end
    
    ipopt_options.ipopt.acceptable_iter = 100;
    ipopt_options.ipopt.tol = 1e-9;
    ipopt_options.ipopt.compl_inf_tol = 1e-6;
    ipopt_options.ipopt.print_level = 0;
    ipopt_options.ipopt.dual_inf_tol = 1e-4;
    ipopt_options.ipopt.constr_viol_tol = 1e-11;
    ipopt_options.ipopt.nlp_scaling_method = 'gradient-based';

    [x, info] = ipopt(x0, f_ipopt, ipopt_options);
    lambda = info.lambda;
    if info.status <= -10
        'Debug';
    end
    
end