function my_checker_min_tr(input, output)

    global accumulator

    for k = 1:numel(accumulator)
        if equal_cells(accumulator{k}.input, input)
            if ~equal_cells(accumulator{k}.output, output)
                warning('cmg:changed', 'Attention');
            end
        end
    end
    accumulator{end+1}.input = input;
    accumulator{end}.output = output;

end