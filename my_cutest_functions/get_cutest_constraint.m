function [c, gc, Hc] = get_cutest_constraint(x, m, factor)
global problem_data_cutest

n_constraints = problem_data_cutest.m;
n_variables = problem_data_cutest.n;

if nargin < 3 || isempty(factor)
    factor = 1;
end

if m <= n_constraints
    if nargout <= 1
        c = cutest_cons(x, m);
    else
        [c, gc] = cutest_cons(x, m);
    end
    if nargout >= 3
        Hc = cutest_ihess(x, m);
    end
elseif m <= n_constraints + n_variables
    % It must be a lower bound
    k = m - n_constraints;
    c = problem_data_cutest.bl(k) - x(k);
    gc = zeros(n_variables, 1);
    gc(k) = -1;
    Hc = zeros(n_variables);    
else
    % Upper bound, then
    k = m - n_constraints - n_variables;
    c = x(k) - problem_data_cutest.bu(k);
    gc = zeros(n_variables, 1);
    gc(k) = 1;
    Hc = zeros(n_variables);    
end

c = factor*c;
if nargout >= 2
    gc = factor*gc;
end
if nargout >= 3
    Hc = factor*Hc;
end
