
size1 = floor(dim/2);
size2 = dim - size1;
M1 = zeros(dim, size1);
M2 = zeros(dim, size2);
for n = 1:size1
    M1(n, n) = 1;
    M2(size1 + n, n) = 1;
end
M2(size1 + size2, size2) = 1;


x1 = M1'*x0;
x2 = M2'*x0;
%%
for k = 1:10
    p1 = @(x) l1_function(f, all_con, mu, M1*x + M2*x2);
    x1 = fmincon(p1, x1,[],[],[],[],[],[], [], fmincon_options);
    p2 = @(x) l1_function(f, all_con, mu, M2*x + M1*x1);
    x2 = fmincon(p2, x2,[],[],[],[],[],[], [], fmincon_options);
end
