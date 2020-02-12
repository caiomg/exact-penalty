function print_result(result, fd)
    if nargin < 2 || isempty(fd)
        fd = 1;
    end
    fprintf(fd, '%s\n', format_result(result));
end
