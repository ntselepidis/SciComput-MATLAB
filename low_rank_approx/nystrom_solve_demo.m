clear; clc; rng(0);

n = 1000;
X = rand(n, 2);
t = ones(n, 1);
m = 100;
sigma = 1e-1;
k = @(x, y) exp(-pdist2(x, y, 'squaredeuclidean'));
alpha = nystrom_solve(X, t, m, sigma, k, true);
