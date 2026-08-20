% Example usage in Octave / MATLAB
%n = 1000;
%A = sprandsym(n, 0.01);   % random sparse symmetric matrix
nx = 32;
n = nx*nx;
A = gallery('poisson', nx);
v1 = randn(n,1);
m = 30;                            % subspace dimension

[V, H] = arnoldi(A, v1, m);

% Small projected matrix
Hm = H(1:m, 1:m);

% Ritz eigenvalues (approximate eigenvalues of A)
ritz = eig(Hm);
