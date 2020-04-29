clear; clc; close all;
[a, b] = deal(-1, 1);
m = 32; % number of subintervals per dimension
f = 5e8; % if f = 0 then Helmholtz -> Poisson
prec = 1; % enable/disable preconditioning
tol = 1e-10;
Nmax = 1000;
restart = [10 20 40 80];
fem_type = 1; % choose 0 for triangles or 1 for squares
if ( fem_type )
    [h, ne, n, coo, con, bounds] = qfem(a, b, m);
    [A, b] = qfem_assemble(coo, con, f);
else
    [h, ne, n, coo, con, bounds] = tfem(a, b, m);
    [A, b] = tfem_assemble(coo, con, f);
end
bounds = [1:(m+1) (m+2):(m+1):((m+1)^2-m)]; % Dirichlet(0) on NE
in = setdiff(1:n, bounds);
xx = zeros(n, 1);
A(bounds, :) = [];
A(:, bounds) = [];
b(bounds) = [];
ex = A \ b; % Direct solution
xx(in) = ex;
xx(bounds) = 0;
figure, trimesh(con, coo(:,1), coo(:,2), xx), % plot solution
title(sprintf('Solution of Helmholtz PDE for f = %e', f));

fprintf(' ------------ Matrix statistics ----------- \n');
fprintf('length(A) = %d\n', length(A));
fprintf('condest(A) = %e\n', condest(A));
eigvals = eig(A); mineig = min(eigvals); maxeig = max(eigvals);
% figure, stem(1:length(A), eigvals), 
% title(['Eigenvalues of A are in D = [' ...
%     num2str(mineig) ',' num2str(maxeig) ']']);

% Testing Krylov Solvers (cg, bicgstab, gmres) and
% Preconditioners (Jacobi, Gauss-Seidel, Symmetric Gauss-Seidel)
D = diag(diag(A));
L = tril(A);

if ( (maxeig*mineig) > 0 ) % check for same sign
    fprintf('\n ------------- cg convergence ------------- \n');
    x = pcg(-A, -b, tol, Nmax); % A is negative definite
    if (prec)
        x = pcg(-A, -b, tol, Nmax, @(y) -D\y);
        x = pcg(-A, -b, tol, Nmax, @(y) -L\y);
        x = pcg(-A, -b, tol, Nmax, @(y) -L\(D*(L'\y)));
    end
end

fprintf('\n ---------- bicgstab convergence ---------- \n');
Nmax2 = ceil(Nmax/2);
x = bicgstab(A, b, tol, Nmax2);
if (prec)
    x = bicgstab(A, b, tol, Nmax2, @(y) D\y);
    x = bicgstab(A, b, tol, Nmax2, @(y) L\y);
    x = bicgstab(A, b, tol, Nmax2, @(y) L\(D*(L'\y)));
end

for i = 1:length(restart)
    fprintf('\n ---------- gmres(%d) convergence --------- \n', restart(i));
    Nmax3 = ceil(Nmax/restart(i));
    x = gmres(A, b, restart(i), tol, Nmax3);
    if (prec)
        x = gmres(A, b, restart(i), tol, Nmax3, @(y) D\y);
        x = gmres(A, b, restart(i), tol, Nmax3, @(y) L\y);
        x = gmres(A, b, restart(i), tol, Nmax3, @(y) L\(D*(L'\y)));
    end
end
