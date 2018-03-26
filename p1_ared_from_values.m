function ared = p1_ared_from_values(fvalues_old, fvalues_new, mu, ind_eactive, multipliers)

if nargin < 4
    ind_eactive = [];
end
if nargin < 5
    multipliers = zeros(size(ind_eactive));
end

n_constraints = size(fvalues_old, 1) - 1;
pos_not_eactive = true(n_constraints, 1);
pos_not_eactive(ind_eactive) = false(size(ind_eactive));

f_change = fvalues_old(1) - fvalues_new(1);
phi_values_old = fvalues_old(2:end);
phi_values_new = fvalues_new(2:end);

m_change = multipliers.*(phi_values_old(ind_eactive) - phi_values_new(ind_eactive));

positive_to_positive = pos_not_eactive & (phi_values_old > 0) & (phi_values_new > 0);
positive_to_negative = pos_not_eactive & (phi_values_old > 0) & (phi_values_new <= 0);
negative_to_positive = pos_not_eactive & (phi_values_old <= 0) & (phi_values_new > 0);
negative_to_negative = pos_not_eactive & (phi_values_old <= 0) & (phi_values_new <= 0);

c_change(positive_to_positive) = phi_values_old(positive_to_positive) - phi_values_new(positive_to_positive);
c_change(positive_to_negative) = phi_values_old(positive_to_negative);
c_change(negative_to_positive) = -phi_values_new(negative_to_positive);
c_change(negative_to_negative) = zeros(sum(negative_to_negative), 1);
sum_c = sum(sort(c_change));


ared = sum(sort([f_change, mu*sum_c, m_change']));

end