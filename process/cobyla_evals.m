function n_evals = cobyla_evals(summary, solutions)

    ncases = length(summary);
    n_evals = nan(ncases, 1);
    max_error = 0;
    max_lgrad = 0;
    for k = 1:ncases
        if max(summary{k}.convals) < 1e-1
            tnorm = norm(summary{k}.lgrad);
            f_error = solutions(k) - summary{k}.fx;
            if tnorm < 1e-1 || -f_error < 1e-4
                n_evals(k) = summary{k}.fcount;
                if max_lgrad < tnorm && tnorm < 1e-1
                    max_lgrad = tnorm;
                end
                if max_error < -f_error && -f_error < 1e-4
                    max_error = f_error;
                end
            end
        end
    end
    max_lgrad
    max_error
end