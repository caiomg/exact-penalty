get_cobyla_results
cobyla_results = verify_results(cobyla_results);
cobyla_evaluations = cobyla_evals(cobyla_results, [selected_problems.solution]);

get_nomad_results
nomad_results = verify_results(nomad_results);
nomad_evaluations = cobyla_evals(nomad_results, [selected_problems.solution]);



cobyla_succeeded = isfinite(cobyla_evaluations);
sum(cobyla_succeeded)
sum(cobyla_evaluations(cobyla_succeeded))