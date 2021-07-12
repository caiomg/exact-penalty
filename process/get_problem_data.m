function [dim, v_bounds, n_constraints, bound_feasible, nl_feasible] = ...
    get_problem_data(name)

    directory = '~/code/my_problems/';

    [~, prob_interface] = setup_cutest_problem(name, directory);
    dim = prob_interface.dim;
    lower_bounds = prob_interface.lb > -1e18;
    upper_bounds = prob_interface.ub < 1e18;
    v_bounds = sum(lower_bounds) + sum(upper_bounds);
    lc_bounds = prob_interface.con_lb > -1e18;
    uc_bounds = prob_interface.con_ub < 1e18;
    n_constraints = sum(lc_bounds) + sum(uc_bounds);
    bound_feasible = 0 == sum(prob_interface.lb > prob_interface.x0 ...
                              | prob_interface.ub < prob_interface.x0);
    nl_feasible = true;
    for k = 1:prob_interface.n_constraints
        ck = prob_interface.evaluate_constraint(prob_interface.x0, k);
        if prob_interface.con_lb(k) > ck ...
                || prob_interface.con_ub(k) < ck
            nl_feasible = false;
        end
    end
end