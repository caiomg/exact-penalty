
dimension = 5;
H = zeros(dimension);
g = zeros(dimension, 1);
c = dimension^2 + dimension + 1;
v = 0;
for k = 1:dimension
    for k2 = k:dimension
        v = v + 1;
        H(k, k2) = v;
        H(k2, k) = v;
    end
end
for k = 1:dimension
    g(k) = H(end) + k;
end

H = 0*H;
f = @(x) quadratic(H, g, c, x);


% Initializing model structure
model.points = dimension*ones(dimension, 1);
model.fvalues = f(model.points(:, 1));
model.radius = 1e-6;
model.basis =  natural_basis(dimension);


options = struct('tol_radius', 1e-5, 'tol_f', 1e-6, 'eps_c', 1e-5, ...
                 'eta_0', 0, 'eta_1', 0.1, 'gamma_inc', 2, 'gamma_dec', ...
                 0.5, 'initial_radius', 1, 'radius_max', 1e3, ...
                 'criticality_mu', 100, 'criticality_beta', 10, ...
                 'criticality_omega', 0.5, 'basis', 'full quadratic', ...
                 'pivot_threshold', 1/6);



% Completing set of interpolation points and calculating polynomial
% model
model = complete_interpolation_set(model, {f}, options);

[c1, g1, H1] = get_model_matrices(model, 0)

g0 = g1 - H1*model.points(:, 1)
c0 = c1 - g1'*model.points(:, 1) + 0.5*(model.points(:, 1)'*H1* ...
                                        model.points(:, 1))