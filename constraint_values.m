function cvalues = constraint_values(constraints, x)

cvalues = evaluate_constraints(constraints, x);
cvalues = [cvalues.c];

end