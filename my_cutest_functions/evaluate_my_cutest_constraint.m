function [c, gc, Hc]  = evaluate_my_cutest_constraint(x, m, factor)

global all_cutest_constraints

if nargin <= 2 || isempty(factor)
    factor = 1;
end

if nargout <= 1
    c = all_cutest_constraints(m).fun(x);
    c = (c - all_cutest_constraints(m).bound)*all_cutest_constraints(m).sign;
elseif nargout <= 2
    [c, gc] = all_cutest_constraints(m).fun(x);
    c = (c - all_cutest_constraints(m).bound)*all_cutest_constraints(m).sign;
    gc = gc*all_cutest_constraints(m).sign;
elseif nargout <= 3
    [c, gc, Hc] = all_cutest_constraints(m).fun(x);
    c = (c - all_cutest_constraints(m).bound)*all_cutest_constraints(m).sign;
    gc = gc*all_cutest_constraints(m).sign;
    Hc = Hc*all_cutest_constraints(m).sign;
end

    
