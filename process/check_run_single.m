function succeeded = check_run_single(results)

    tol_c = 1e-2;
    tol_f = 1e-1;

    n_problems = size(results, 1);
    succeeded = false(n_problems, 1);
    for k = 1:n_problems
        if ~isempty(results(k, 1).nphi) && results(k, 1).nphi < tol_c
            succeeded(k, 1) = true;
        end
    end

end
