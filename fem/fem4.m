clear; clc; close all;
[a, b] = deal(-1, 1);
m = 64; % number of subintervals per dimension
[h, ne, n, coo, con, bounds] = qfem(a, b, m);
bounds = [1:(m+1) (m+2):(m+1):((m+1)^2-m)]; % Dirichlet (0) on SW bounds, Neumann (0) elsewhere
f = 5e8;
c = 3e8;
k = (2*pi*f)/c;
[A, b] = qfem_assemble(coo, con, f);
A(bounds, :) = [];
A(:, bounds) = [];
b(bounds) = [];
x = A \ b;
in = setdiff(1:n, bounds);
out = bounds;
xx = zeros(n, 1);
xx(in) = x;
xx(out) = 0;
figure, trimesh(con, coo(:, 1), coo(:, 2), xx);
