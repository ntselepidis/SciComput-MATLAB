clear; clc; close all;

n = 10000; % Number of total samples
m = 100;   % Number of samples to be used in approximation
kernel = @(x, y) exp(-pdist2(x, y, 'squaredeuclidean')); % Kernel function
sigma = 0.1; % Jitter factor

% Generate dataset and setup solution of linear system
X = rand(n, 2);
K = kernel(X, X);
I = speye(n);
a = (1:n)' / n;
t = (K + sigma * I) * a;

% Generate random permutation to shuffle dataset
perm = randperm(n);

% Perform Nystrom preprocessing
[Lambda, U] = nystrom_prep(X(perm, :), kernel, m);

% Compute Nystrom approximaton
ahat = nystrom_solve(t(perm), sigma, Lambda, U);

% Back-permute solution
ahat(perm) = ahat;

% Compute error of approximation
relerr = norm(a - ahat) / norm(a)
relres = norm(t - (K + sigma * I) * ahat) / norm(t)

% Plot result
figure, hold on, grid,
plot(a, 'b-');
plot(ahat, 'ro');
