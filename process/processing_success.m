
for k = 1:50
    if ~strcmp(summary(k).name, all_feasibility_info{k}.name) ...
            || ~strcmp(selected_problems(k).name, summary(k).name)
        error();
    else
        if isfield(all_feasibility_info{k}, 'fx')
            decrease_obtained = all_feasibility_info{k}.fx - summary(k).fx;
            maximum_decrease = all_feasibility_info{k}.fx - selected_problems(k).solution;
            decrease_ratio(k, 1) = decrease_obtained/maximum_decrease;
        else
            decrease_ratio(k, 1) = -inf;
        end
        grad_norm(k, 1) = norm(summary(k).lgrad);
        error_obj(k, 1) = summary(k).error_obj;
        viol(k, 1) = summary(k).nphi;
    end
end

successes = (viol < 1e-5 & (decrease_ratio > 1-1e-5 | grad_norm < 1e-5 | -error_obj < 1e-6))
failures = find(~successes)
sum(successes)

summary_corrected = summary;
for k = 1:numel(failures)
    summary_corrected(failures(k)).fcount = nan;
    summary_corrected(failures(k)).best_fcount = nan;
end