clear; close all;

ndoms = 32;  % Number of subdomains
maxit = 500; % Maximum number of iterations for Krylov method
tol = 1e-10; % Prescribed tolerance

%nx = 10;
nx = 32;
A = gallery('poisson', nx);
%A = -nine_point(nx+1);
%A = poisson3D(nx);
%ex = (1:length(A))'/length(A);
%b = A*ex;
b = ones(length(A), 1);

% Partitioning
G = 0.5*(abs(A)+abs(A)'); % undirected graph
map = PartSparseMat(G,ndoms);
ndoms = max(map);

% Classify subdomain components into inner and outer points
[In, Out, All, blk] = inner_outer(G, map);

% Block Jacobi preprocessing
BJ.blk = blk;
[BJ.L, BJ.U, BJ.P, BJ.Q] = factorize(A, All);

% Permute matrix A in a block diagonal form
perm = horzcat( All{:} )';
AP = A(perm, perm);
bP = b(perm);

% Nicolaides coarse space for block Jacobi
V = aggregate(AP, blk);
VAV = V*AP*V';

%VAV2 = zeros(ndoms, ndoms);
%for i = 1:ndoms
%    ii = blk(i):blk(i+1)-1;
%    for j = 1:ndoms
%        jj = blk(j):blk(j+1)-1;
%        VAV2(i,j) = sum(AP(ii,jj), 'all') / (sqrt(length(ii)) * sqrt(length(jj)));
%    end
%    Vb(i) = sum(bP(ii), 'all') / sqrt(length(ii));
%end
%max(abs(sparse(VAV2) - VAV), [], 'all')
%max(abs(Vb - V*b), [], 'all')
%
%(V*V')

% Preconditioned conjugate gradient on Ax = b, with preconditioner M
%x = pcg( AP, bP, tol, maxit );                                          % No preconditioner
%x = pcg( AP, bP, tol, maxit, @(y) blkjac(y,BJ) );                       % Block Jacobi
%x = pcg( AP, bP, tol, maxit, @(y) blkjac(y,BJ) + V'*(VAV\(V*y)) );      % Block Jacobi + Nico
%x = pcg( AP, bP, tol, maxit, @(y) blkjac(y,BJ) + 0.1*V'*(VAV\(V*y)) );  % Block Jacobi + damping * Nico

% Preconditioned Richardson iteration (gradient descent)
% x = x + M*(b - A*x), M is the preconditioner
% x = (I-M*A)*x + M*b, alternative formulation
% here it is written in 3-steps as
% r = b - A*x
% e = M*r
% x = x + e

%figure,
dmp = 0.1;
% poisson2D(32)
%omega = 1; % BJ
omega = 0.6376; % BJ+CGC
% 9-point(32)
%omega = 1.0422; % BJ
%omega = 0.6534; % BJ+CGC
% poisson3D(10)
%omega = 1.0880; % BJ
%omega = 0.7036; % BJ+CGC

%omega = 2./(maxev + minev);

% add-then-invert (woodbury)
Vt_vec = full(sum(V', 2));
Vt_pat = spones(V');
VAV_zero_diag = VAV - diag(diag(VAV));
VAVc = speye(ndoms) + VAV_zero_diag * ( V*( blkjac(Vt_vec, BJ).*Vt_pat ) );

xP = zeros(length(A), 1);
disp(['omega = ' num2str(omega)])
for i = 1 : 500
    rP = bP - AP*xP;
    omega = (rP'*rP) / (rP'*AP*rP);
    relres = norm(rP) / norm(bP);
    disp([i omega relres])
    if (relres <= 1e-3)
        disp(['Converged in ' num2str(i) ' iterations.'])
        break;
    end
    %eP = rP;
    %eP = blkjac(rP, BJ);
    %eP = V'*(VAV\(V*rP));
    eP = blkjac(rP, BJ) + (V'*(VAV\(V*rP)));
    %eP = blkjac(rP, BJ) +  dmp * (V'*(VAV\(V*rP)));
    %y = blkjac(rP, BJ);
    %eP = y - blkjac(V'*(VAVc \ (VAV_zero_diag * (V*y))), BJ);
    xP = xP + omega * eP;
    %x(perm) = xP;
    %getframe; mesh(reshape(x, nx, nx))
end

return
% Don't continue analysis if matrix is too large
assert(length(A) <= 1024);

% Compute preconditioned matrix MA for each case
A = full(A);
MA_bj = zeros(size(A));
MA_cgc = zeros(size(A));
for i = 1 : length(A)
    rP = A(perm, i);
    MA_bj(perm, i) = blkjac(rP, BJ);     % Block Jacobi
    MA_cgc(perm, i) = (V'*(VAV\(V*rP))); % CGC
end
MA_bj_cgc  = MA_bj + MA_cgc;             % Block Jacobi + CGC
MA_bj_dcgc = MA_bj + dmp * MA_cgc;       % Block Jacobi + damped CGC

% Compute eigenvalues of preconditioned matrix MA for each case
v         = eig(A);
v_bj      = eig(MA_bj);
v_cgc     = eig(MA_cgc);
v_bj_cgc  = eig(MA_bj_cgc);
v_bj_dcgc = eig(MA_bj_dcgc);

% Print some statistics
disp('Condition number of MA')
kappa = [cond(A), cond(MA_bj), cond(MA_cgc), cond(MA_bj_cgc), cond(MA_bj_dcgc)]

disp('Max, Min and, Max/Min Ratio of Eigenvalues of MA')
maxev = [max(abs(v)), max(abs(v_bj)), max(abs(v_cgc)), max(abs(v_bj_cgc)), max(abs(v_bj_dcgc))]
minev = [min(abs(v)), min(abs(v_bj)), min(abs(v_cgc)), min(abs(v_bj_cgc)), min(abs(v_bj_dcgc))]
ratio = maxev ./ minev

disp('Optimal omega')
omega = 2./(maxev + minev)

I = speye(size(A));
disp('norm(I-MA)')
nrm   = [norm(I-A), norm(I-MA_bj), norm(I-MA_cgc), norm(I-MA_bj_cgc), norm(I-MA_bj_dcgc)]

disp('Max Eigenvalue of (I-MA)')
maxevIminusMA = [ ...
    max(abs(eig(I-A))), ...
    max(abs(eig(I-MA_bj))), ...
    max(abs(eig(I-MA_cgc))), ...
    max(abs(eig(I-MA_bj_cgc))), ...
    max(abs(eig(I-MA_bj_dcgc))) ...
    ]

% return

figure,

subplot(5, 1, 1),
plot(real(v), imag(v), 'b*'),
title(['eig(A), cond(A)  = ' ...
    num2str(kappa(1)) ' \geq \lambda_{max}/\lambda_{min} = ' num2str(ratio(1))]);

subplot(5, 1, 2),
plot(real(v_bj), imag(v_bj), 'ro'),
title(['eig(MA), M = BJ, cond(MA) = ' ...
    num2str(kappa(2)) ' \geq \lambda_{max}/\lambda_{min} = ' num2str(ratio(2))]);

subplot(5, 1, 3),
plot(real(v_cgc), imag(v_cgc), 'go'),
title(['eig(MA), M = CGC, cond(MA) = ' ...
    num2str(kappa(3)) ' \geq \lambda_{max}/\lambda_{min} = ' num2str(ratio(3))]);

subplot(5, 1, 4),
plot(real(v_bj_cgc), imag(v_bj_cgc), 'cx'),
title(['eig(MA), M = BJ+CGC, cond(MA) = ' ...
    num2str(kappa(4)) ' \geq \lambda_{max}/\lambda_{min} = ' num2str(ratio(4))]);

subplot(5, 1, 5),
plot(real(v_bj_dcgc), imag(v_bj_dcgc), 'kx'),
title(['eig(MA), M = BJ+dmp*CGC, cond(MA) = ' ...
    num2str(kappa(5)) ' \geq \lambda_{max}/\lambda_{min} = ' num2str(ratio(5))]);

% return

% Visualization
% Exact Solution
figure, mesh(reshape(A\b, nx, nx));
% Block Jacobi
xP = blkjac(bP, BJ);
x(perm) = xP;
figure, mesh(reshape(x, nx, nx));
% Coarse Grid
xP = V'*(VAV\(V*bP));
x(perm) = xP;
figure, mesh(reshape(x, nx, nx));
