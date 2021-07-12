function result = strincell(c, s)

    result = false;
    for k = 1:numel(c)
        if strcmp(c{k}, s)
            result = true;
            break
        end
    end
end