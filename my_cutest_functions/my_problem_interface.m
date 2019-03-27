classdef my_problem_interface
    
    properties
        directory
        n_constraints
    end
    
    methods
        function [f, g, H] = evaluate_objective(self, x)
            old_dir = cd(self.directory);
            cleanup_obj = onCleanup(@() cd(old_dir));
            if nargout == 1
                f = get_cutest_objective(x);
            elseif nargout == 2
                [f, g] = get_cutest_objective(x);
            else
                [f, g, H] = get_cutest_objective(x);
            end
        end
        function [f, g, H] = evaluate_constraint(self, x, m)
            assert(m <= self.n_constraints);
            old_dir = cd(self.directory);
            cleanup_obj = onCleanup(@() cd(old_dir));
            if nargout == 1
                f = get_cutest_constraint(x, m);
            elseif nargout == 2
                [f, g] = get_cutest_constraint(x, m);
            else
                [f, g, H] = get_cutest_constraint(x, m);
            end
        end
    end
    
end
