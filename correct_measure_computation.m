function [measure, d, theta] = correct_measure_computation(pg, G, mu, lb, ub, dt)
% CORRECT_MEASURE_COMPUTATION - 
%   
    
    dim = numel(pg);
    n_aconstraints = size(G, 2);
    assert (dim == numel(lb) && dim == numel(ub));
    
    if isempty(dt)
        measure = [];
        d = [];
        theta = [];
    else
        assert (numel(dt) == dim + n_aconstraints);
        d_ot = dt(1:dim);
        theta_ot = dt(dim+1:dim+n_aconstraints);
        
        d = project_to_bounds(d_ot, lb, ub);
        if n_aconstraints > 0
            theta = max(0, G'*d);
        else
            theta = 0;
        end
        
        measure = -(pg'*d + mu*norm(theta, 1));
        
        error_d = d_ot - d;
        error_theta = theta_ot - theta;        
    end
    
end
