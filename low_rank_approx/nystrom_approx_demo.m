clear; clc; rng(0);

n = 10000;
X = rand(n, 2);
t = (1:n)' / n;
m = 100;
sigma = 1e-1;
compute_error = true;

% Kernel function
kernel = @(x, y) exp(-pdist2(x, y, 'squaredeuclidean'));

% Shuffle dataset
perm = randperm(n);
Xperm = X(perm, :);
tperm = t(perm);

% Perform Nystrom preprocessing
[Lambda, U] = nystrom_prep(Xperm, kernel, m);

% Compute Nystrom approximaton
alpha = nystrom_solve(tperm, sigma, Lambda, U);

% Back-permute solution
alpha(perm) = alpha;

if compute_error
    
    % Compute full kernel matrix
    Knn = kernel(X, X);

    % Solve full linear system
    alpha_true = (Knn + sigma * speye(size(X, 1))) \ t;

    % Compute error of approximation
    disp(norm(alpha - alpha_true) / norm(alpha_true));

end
