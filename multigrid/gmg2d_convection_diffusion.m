clear; clc; close all;
lev = 6; % h = 1/(2^lev)
mu = 1; % choose 1 for V-cycle or 2 for W-cycle
str = 'VW';
omega = 4/5;
% omega = 1;
tol = 1e-10;
Nmax = 50;
n1 = 2; % Number of pre-smoothing iterations
n2 = 2; % Number of post-smoothing iterations
epsilon = 1; % Diffusion Coefficient
alpha = 25; % Convection Coefficient
A = cell(lev, 1); % Coefficient Matrix
T = cell(lev, 1); % Diffusion Matrix
K = cell(lev, 1); % Convection Matrix
M = cell(lev, 1); % Smoother
b = cell(lev, 1); % Rhs vector (source)
x = cell(lev, 1); % Solution vector
P = cell(lev-1, 1); % Prolongation Matrix
R = cell(lev-1, 1); % Restriction Matrix
for i = 1:lev
    n = 2^(lev-i+1)-1;
    T{i} = epsilon*gallery('tridiag', n)*(n+1)*(n+1);
%     K{i} = alpha*spdiags([-ones(n, 1) ones(n, 1)], [-1 1], n, n)*((n+1)/2);
    K{i} = alpha*spdiags([-ones(n, 1) ones(n, 1)], [-1 0], n, n)*(n+1);
    A{i} = kron(T{i}, speye(n))+kron(speye(n), T{i}+K{i});
    M{i} = diag(diag(A{i})); % Jacobi
%     M{i} = tril(A{i}); % Gauss-Seidel
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
b{1} = ones(length(A{1}), 1);
disp(['2D Convection-Diffusion Equation : ' num2str(length(A{1}))]);
disp('Relative residual');disp(' ');
for i = 1:Nmax
    x = mgmu3(A, b, x, P, R, M, n1, n2, 1, lev, omega, mu); % Implicit Smoothing
    normr = norm(b{1}-A{1}*x{1}) / norm(b{1});
    disp(normr);
    if (normr < tol )
        disp([str(mu) '-cycle multigrid converged in ' ...
            num2str(i) ' iterations.']);
        break;
    end 
end
nx = sqrt(length(x{1}));
row = A{1}(floor((nx^2)/2), :);
row = full(row(row ~= 0)); disp(' '); disp(row);
ddom = abs(row(3))-sum(abs(row([1 2 4 5]))) % diagonal dominance
condest(A{1})
figure,  mesh(reshape(x{1}, nx, nx)), 
title(['Multigrid solution for epsilon = ' ...
    num2str(epsilon) ' and alpha = ' num2str(alpha)]);
figure,  mesh(reshape(A{1}\b{1}, nx, nx)), 
title(['Direct solution for epsilon = ' ...
    num2str(epsilon) ' and alpha = ' num2str(alpha)]);
