function print_oogp_table(results)
    
    print_row(results);
    
end



function print_row(results)
    fprintf(1, '%g, %g, %g\n', results.qo, results.qg, results.zco2);
end
