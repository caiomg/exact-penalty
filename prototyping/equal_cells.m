function result = equal_cells(a, b)

    if numel(a) ~= numel(b)
        result = false;
    else
        result = true;
        for k = 1:numel(a)
            if ~isempty(find(a{k} ~= b{k}, 1))
                result = false;
            end
        end
    end
end
