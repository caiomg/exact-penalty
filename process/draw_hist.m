function draw_hist(summary)

    ncases = length(summary);
    lag_norms = [];
    for k = 1:ncases
        if max(summary{k}.convals) < 1e-1
            tnorm = norm(summary{k}.lgrad);
            if tnorm < 1e-1
                lag_norms(end+1) = tnorm;
            end
        end
    end
    histogram(lag_norms)

end