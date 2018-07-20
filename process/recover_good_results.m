
clear new_good_results
rfound = 0;
for k = 1:length(all_results)
    this_results = all_results{k};
    for m = 1:length(this_results)
        if ~isempty(this_results(m).fx)
            if this_results(m).kkt || (this_results(m).nphi < 1e-6 && this_results(m).error_obj < 1e-6)
                rfound = rfound+1;
                new_good_results(rfound, 1) = this_results(m);
            end
        end
    end
end

            
    
    