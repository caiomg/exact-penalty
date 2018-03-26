function p = l1_penalty_from_values(fvalues, mu, ind_eactive, multipliers)

if nargin < 4
    multipliers = [];
end
if nargin < 3
    ind_eactive = [];
end

n_constraints = size(fvalues, 1) - 1;
ind_not_eactive = (1:n_constraints)';
ind_not_eactive(ind_eactive) = [];

violations_val = max(0, fvalues(ind_not_eactive + 1));
activities_val = fvalues(ind_eactive + 1);

p = fvalues(1) + mu*sum(violations_val) + multipliers'*activities_val;

end