clear; clc; close all;
lev = 7;
mu = 1;
str = 'VW';
omega = 1;% Needs DOUR relaxation to compute omega
tol = 1e-10;
Nmax = 50;
n1 = 2;
n2 = 2;
A = cell(lev, 1);
D = cell(lev, 1);
b = cell(lev, 1);
x = cell(lev, 1);
P = cell(lev-1, 1);
R = cell(lev-1, 1);
nx = 50;
A{1} = gallery('poisson', nx);
n = length(A{1});
b{1} = A{1}*(1:n)';
x{1} = zeros(n, 1);
lfill = 1;
D{1} = make_spai(A{1}, lfill);
for i = 2:lev
    R{i-1} = coarsening(A{i-1});
    P{i-1} = R{i-1}';
    A{i} = R{i-1}*A{i-1}*P{i-1};
    D{i} = make_spai(A{i}, lfill);
    x{i} = zeros(length(A{i}), 1);
    b{i} = zeros(length(A{i}), 1);
end
disp(['2D Poisson : ' num2str(length(A{1}))]);
disp('Relative residual');disp(' ');
for i = 1:Nmax
    x = mgmu2(A, b, x, P, R, D, n1, n2, 1, lev, omega, mu);
    normr = norm(b{1}-A{1}*x{1}) / norm(b{1});
    disp(normr);
    if ( normr < tol )
        disp([str(mu) '-cycle multigrid converged in ' num2str(i) ' iterations.']);
        break;
    end 
end
