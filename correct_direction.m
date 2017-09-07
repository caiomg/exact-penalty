function h = correct_direction(h, A)

if size(h, 1) ~= size(A, 1)
    error();
elseif size(h, 2) ~= 1
    error();
end

for column = A
    h = h - column*(h'*column)/(column'*column);
end

inexact_coordinates = find(A'*h < 0);
if ~isempty(inexact_coordinates)
    for n = inexact_coordinates'
        h = h - 2*A(:, n)*(h'*A(:, n))/(A(:, n)'*A(:, n));
    end
end


end