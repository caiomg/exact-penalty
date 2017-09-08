function [ind_eactive, ind_eviolated] = ...
                     identify_new_constraints(current_constraints, epsilon, ...
                                              ind_eactive)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% ind_eactive = [];
ind_eviolated = zeros(0, 1);
n_constraints = size(current_constraints, 1);
for n = 1:n_constraints
    if isempty(find(ind_eactive == n, 1))
        if abs(current_constraints(n).c) <= epsilon
            ind_eactive(end+1, 1) = n;
        % Should this be here?
        % Maybe one constraint could be active and included as violated
        elseif current_constraints(n).c < epsilon
            ind_eviolated(end+1, 1) = n;
        end
    end
end

end

