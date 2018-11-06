function good_results = check_good_results(results)

    good_results = {};
    all_errors= [];
    for k = 1:length(results)
        if results(k).kkt || ...
            (~isempty(results(k).nphi) && results(k).nphi < 1e-6 && ...
             (...%-results(k).error_obj < 1e-5 || ...
                -results(k).error_obj/abs(results(k).fx + results(k).error_obj) < 1e-6))
            good_results{end+1} = results(k);
        end
    end
        
end