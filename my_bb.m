function [d, upper_included, fval] = my_bb(p_grad, mu, A, lb, ub, fixed_constraints, best)
    
% Constraints enforced
    fixed_lower = fixed_constraints < 0;
    A_lower = A(fixed_lower);
    fixed_upper = fixed_constraints > 0;
    A_upper = A(fixed_upper);
    Aineq = [A_lower;
             -A_upper];
    bineq = zeros(size(Aineq, 1));
    
    f = p_grad + mu*sum(A_upper)';
    
    d = solve_linear_problem(f, Aineq, bineq, [], [], lb, ub);
    
    if f'*d > best
        fval = f'*d; % bound
        upper_included = fixed_upper;
    else
        upper = A*d > 0;
        new_upper = upper & ~fixed_upper;
        next = find(new_upper, 1);
        if ~isempty(next)
            fixed_constraits_a = fixed_constraints;
            fixed_constraints_a(next) = -1;
            [d_a, upper_icluded_a, fval_a] = ...
                my_bb(p_grad, mu, A, lb, ub, fixed_constraints_a, best);
            if fval_a < best
                best = fval_a;
                d = d_a;
                upper_included = upper_included_a;
                fval = fval_a;
            end
            fixed_constraits_b = fixed_constraints;
            fixed_constraints_b(next) = 1;
            [d_b, upper_icluded_b, fval_b] = ...
                my_bb(p_grad, mu, A, lb, ub, fixed_constraints_b, best);
            if fval_b < best
                d = d_b;
                upper_included = upper_included_b;
                fval = fval_b;
            end
        else
            % optimal
            fval = f'*d;
            upper_included = upper;
        end
    end
end