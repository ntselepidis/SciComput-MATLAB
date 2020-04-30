clear; clc; close all;
lev = 7;
mu = 1;
str = 'VW';
omega = 4/5;%2/3;
tol = 1e-10;
Nmax = 50;
n1 = 2;
n2 = 2;
A = cell(lev, 1);
T = cell(lev, 1);
D = cell(lev, 1);
b = cell(lev, 1);
x = cell(lev, 1);
P = cell(lev-1, 1);
R = cell(lev-1, 1);
for i = 1:lev
    n = 2^(lev-i+1)-1;
    T{i} = gallery('tridiag', n)*(n+1)*(n+1);
    A{i} = kron(T{i}, speye(n))+kron(speye(n), T{i});
    D{i} = diag(A{i});
    x{i} = zeros(n*n, 1);
    b{i} = zeros(n*n, 1);
end
for i = 1:lev-1
    r = 2^(lev-i+1)-1;
    c = 2^(lev-i)-1;
    P{i} = spalloc(r, c, 3*c);
    for j = 1:c
        P{i}((1:3)+2*(j-1), j) = [0.5 1 0.5]';
    end
    R{i} = 0.5*P{i}';
    P{i} = kron(P{i}, P{i});
    R{i} = kron(R{i}, R{i});
end
b{1} = A{1}*(1:((2^lev-1)^2))';
disp(['2D Poisson : ' num2str(length(A{1}))]);
disp('Relative residual');disp(' ');
for i = 1:Nmax
    x = mgmu(A, b, x, P, R, D, n1, n2, 1, lev, omega, mu);
    normr = norm(b{1}-A{1}*x{1}) / norm(b{1});
    disp(normr);
    if ( normr < tol )
        disp([str(mu) '-cycle multigrid converged in ' num2str(i) ' iterations.']);
        break;
    end 
end
