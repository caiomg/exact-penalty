function faz_grafico(f, points)

n_points = size(points, 2);
y = zeros(n_points, 1);
for k = 1:n_points
    y(k) = f(points(:, k));
end
plot(y);
end