function model = discard_points(model, discard_radius)
% DISCARD_POINTS removes from set of interpolation points that are
% far from center. Does not update polynomial model

center = model.points(:, 1);
points = model.points;
[dim, n_points] = size(points);

tol = dim*10*eps(max(1, discard_radius));

keep = true(1, n_points); % points to keep
for k = 2:n_points
    if (norm(points(:, k) - center, 2) - discard_radius > tol)
        keep(k) = false;
    end
end

% n_notkeep = sum(~keep);
% model.cached_points(:, end+1:end+n_notkeep) = points(:, ~keep);
% model.cached_fvalues(:, end+1:end+n_notkeep) = model.fvalues(:, ~ ...
%                                                   keep);

model.points = points(:, keep);
model.fvalues = model.fvalues(:, keep);


end