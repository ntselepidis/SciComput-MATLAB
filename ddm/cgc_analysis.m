clear; close all;

rand('seed', 0);

ndoms = 32;  % Number of subdomains
maxit = 500; % Maximum number of iterations for Krylov method
tol = 1e-10; % Prescribed tolerance

%nx = 10;
nx = 32;
h = 1 / (nx+1);
A = gallery('poisson', nx) / h^2;
%A = -nine_point(nx+1) / (6*(h^2));
%A = poisson3D(nx) / h^3;
%ex = (1:length(A))'/length(A);
%b = A*ex;
b = ones(length(A), 1);

% Partitioning
G = 0.5*(abs(A)+abs(A)'); % undirected graph
%map = PartSparseMat(G,ndoms);
[x_, y_]     = meshgrid(h:h:1-h); coo = [x_(:) y_(:)];
%[x_, y_, z_] = meshgrid(h:h:1-h); coo = [x_(:) y_(:) z_(:)];
map = kmeans(coo, ndoms);
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

% Eigenvectors of diagonal blocks of matrix AP
%n = length(AP);
%nbasis = 1;
%V = spalloc(nbasis*ndoms, n, nbasis*n);
%for i = 1:ndoms
%    idx = blk(i):blk(i+1)-1;
%    if nbasis > 1
%        [eigvecs, eigvals] = eigs(AP(idx,idx), nbasis-1, 'sm');
%    end
%    V(nbasis*(i-1)+1, idx) = 1;% / length(blk(i):blk(i+1)-1);
%    for j = 2:nbasis
%        V(nbasis*(i-1)+j,idx) = eigvecs(:,j-1);
%    end
%end

% Eigenvectors of full matrix AP
%n = length(AP);
%nbasis = 1;
%[eigvec, eigval] = eigs(AP, nbasis, 'sm');
%V = spalloc(nbasis*ndoms, n, nbasis*n);
%for i = 1:ndoms
%    idx = blk(i):blk(i+1)-1;
%    for j = 1:nbasis
%        V(nbasis*(i-1)+j,idx) = eigvec(idx,j);
%    end
%end

VAV = V*AP*V';

% SHEM Coarse Space for block Jacobi
%nbasis = 1;
%[V, basis] = shem(A, In, Out, All, nbasis);
%V = V(:, perm);
%VAV = V*AP*V';

%VAV2 = zeros(ndoms, ndoms);
%for i = 1:ndoms
%    ii = blk(i):blk(i+1)-1;
%    for j = 1:ndoms
%        jj = blk(j):blk(j+1)-1;
%        VAV2(i,j) = sum(sum(AP(ii,jj))) / (sqrt(length(ii)) * sqrt(length(jj)));
%    end
%    Vb(i) = sum(bP(ii)) / sqrt(length(ii));
%end
%max(max(abs(sparse(VAV2) - VAV)))
%max(max(abs(Vb - V*b)))
%
%(V*V')

% deflation
P  = @(x) x - AP*(V'*(VAV\(V*x)));
Pt = @(x) x - V'*(VAV\(V*(AP*x)));

% Preconditioned richardson iteration on Ax = b, with preconditioner M
disp('')
disp('Preconditioned richardson iteration')
disp('')
maxit = 100;
[x, val0] = richardson( @(y) AP*y, bP, tol, maxit, @(y) y );                                                       % No preconditioner
[x, val1] = richardson( @(y) AP*y, bP, tol, maxit, @(y) blkjac(y, BJ) );                                           % Block Jacobi
[x, val2] = richardson( @(y) AP*y, bP, tol, maxit, @(y) blkjac(y, BJ) + V'*(VAV\(V*y)) );                          % Block Jacobi + CGC (additive)
[x, val3] = richardson( @(y) AP*y, bP, tol, maxit, @(y) blkjac(P(y), BJ) + V'*(VAV\(V*y)) );                       % Block Jacobi + CGC (mult, CGC 1st)
[x, val4] = richardson( @(y) AP*y, bP, tol, maxit, @(y) Pt(blkjac(y, BJ)) + V'*(VAV\(V*y)), [], V'*(VAV\(V*bP)) ); % Block Jacobi + CGC (mult, CGC 2nd)
[x, val5] = richardson( @(y) AP*y, bP, tol, maxit, @(y) Pt(blkjac(P(y), BJ)) + V'*(VAV\(V*y)) );                   % BNN (single-step deflation)
[x, val6] = richardson( @(y) AP*y, bP, tol, maxit, @(y) Pt(blkjac(P(y), BJ)), [], V'*(VAV\(V*bP)) );               % R-BNN1
[x, val7] = richardson( @(y) AP*y, bP, tol, maxit, @(y) Pt(blkjac(y, BJ)), [], V'*(VAV\(V*bP)) );                  % R-BNN2
[x, val8] = richardson( @(y) P(AP*y), P(bP), tol, maxit, @(y) blkjac(y, BJ), [] );                                 % Block Jacobi + deflation
figure, grid, hold,
plot(val0, 'linewidth', 2, 'k')
plot(val1, 'linewidth', 2, 'b')
plot(val2, 'linewidth', 2, 'r')
%plot(val3, 'linewidth', 2, 'm')
%plot(val4, 'linewidth', 2, 'c')
plot(val5, 'linewidth', 2, 'g')
%plot(val6, 'linewidth', 2, 'gx')
%plot(val7, 'linewidth', 2, 'co')
%plot(val8, 'linewidth', 2, 'k--')
return

% Preconditioned conjugate gradient on Ax = b, with preconditioner M
%disp('')
%disp('Preconditioned conjugate gradient')
%disp('')
%[x, ~, ~, ~, val0] = pcg( AP, bP, tol, maxit );                                                               % No preconditioner
%[x, ~, ~, ~, val1] = pcg( AP, bP, tol, maxit, @(y) blkjac(y, BJ) );                                           % Block Jacobi
%[x, ~, ~, ~, val2] = pcg( AP, bP, tol, maxit, @(y) blkjac(y, BJ) + V'*(VAV\(V*y)) );                          % Block Jacobi + CGC (additive)
%[x, ~, ~, ~, val3] = pcg( AP, bP, tol, maxit, @(y) blkjac(P(y), BJ) + V'*(VAV\(V*y)) );                       % Block Jacobi + CGC (mult, CGC 1st)
%[x, ~, ~, ~, val4] = pcg( AP, bP, tol, maxit, @(y) Pt(blkjac(y, BJ)) + V'*(VAV\(V*y)), [], V'*(VAV\(V*bP)) ); % Block Jacobi + CGC (mult, CGC 2nd)
%[x, ~, ~, ~, val5] = pcg( AP, bP, tol, maxit, @(y) Pt(blkjac(P(y), BJ)) + V'*(VAV\(V*y)) );                   % BNN (single-step deflation)
%[x, ~, ~, ~, val6] = pcg( AP, bP, tol, maxit, @(y) Pt(blkjac(P(y), BJ)), [], V'*(VAV\(V*bP)) );               % R-BNN1
%[x, ~, ~, ~, val7] = pcg( AP, bP, tol, maxit, @(y) Pt(blkjac(y, BJ)), [], V'*(VAV\(V*bP)) );                  % R-BNN2
%[x, ~, ~, ~, val8] = pcg( @(x) P(AP*x), P(bP), tol, maxit, @(y) blkjac(y, BJ) );                              % Block Jacobi + deflation
%figure, grid, hold,
%plot(log2(val0), 'linewidth', 2, 'k')
%plot(log2(val1), 'linewidth', 2, 'b')
%plot(log2(val2), 'linewidth', 2, 'r')
%plot(log2(val3), 'linewidth', 2, 'm')
%plot(log2(val4), 'linewidth', 2, 'c')
%plot(log2(val5), 'linewidth', 2, 'g')
%plot(log2(val6), 'linewidth', 2, 'gx')
%plot(log2(val7), 'linewidth', 2, 'co')
%plot(log2(val8), 'linewidth', 2, 'k--')
%return

%x(perm) = V'*(VAV\(V*bP)) + Pt(x);
%mesh(reshape(x, nx, nx))
%
%return

% Preconditioned Richardson iteration (gradient descent)
% x = x + M*(b - A*x), M is the preconditioner
% x = (I-M*A)*x + M*b, alternative formulation
% here it is written in 3-steps as
% r = b - A*x
% rhat = M*r
% x = x + omega*rhat

%figure,
dmp = 2;
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
%VAVc = speye(ndoms) + VAV_zero_diag * ( V*( blkjac(Vt_vec, BJ).*Vt_pat ) );
%VAVc = speye(ndoms) + VAV_zero_diag * ( V*sparse( blkjac(Vt_vec, BJ).*full(Vt_pat) ) );

%prec = @(rP) rP;                                        % No preconditioner
%prec = @(rP) blkjac(rP, BJ);                            % Block Jacobi
%prec = @(rP) V'*(VAV\(V*rP));                           % CGC
prec = @(rP) blkjac(rP, BJ) + (V'*(VAV\(V*rP)));        % Block Jacobi + CGC (additive)
%prec = @(rP) blkjac(P(rP), BJ) + (V'*(VAV\(V*rP)));     % Block Jacobi + CGC (mult, CGC 1st)
%prec = @(rP) Pt(blkjac(rP, BJ)) + (V'*(VAV\(V*rP)));    % Block Jacobi + CGC (mult, CGC 2nd)
%prec = @(rP) Pt(blkjac(P(rP), BJ)) + (V'*(VAV\(V*rP))); % BNN (single-step deflation)
%prec_helper = @(rhatP_) rhatP_ - blkjac(V'*(VAVc \ (VAV_zero_diag * (V*rhatP_))), BJ);
%prec = @(rP) prec_helper(blkjac(rP, BJ));               % Add-then-invert
%mesh = @(x) surf(x);

xP = zeros(length(A), 1);
disp(['omega = ' num2str(omega)])
for i = 1 : 500
    rP = bP - AP*xP;
    %rP = P(rP); % deflation
    relres = norm(rP) / norm(bP);
    disp([i omega relres])
    if (relres <= 1e-3)
        disp(['Converged in ' num2str(i) ' iterations.'])
        break;
    end
    % -------------
    % Batch update
    % -------------
    rhatP = prec(rP);
    %omega = (rhatP'*rhatP) / (rhatP'*(AP*prec(rhatP)));
    omega = (rP'*rhatP) / (rhatP'*AP*rhatP);
    xP = xP + omega * rhatP;
    % -------------
    % Online update
    % -------------
    %omega = 1.0;
    %xP = xP + omega * V'*(VAV\(V*rP));    % Coarse-space correction
    %rP = bP - AP*xP;
    %xP = xP + omega * blkjac(rP, BJ);     % Block Jacobi smoothing
    %xP = blkgs(xP, BJ, @(xP) bP - AP*xP); % Block GS smoothing
    %rP = bP - AP*xP;
    %xP = xP + omega * V'*(VAV\(V*rP));    % Coarse-space correction
    %rP = bP - AP*xP;
    %xP = xP + omega * blkjac(rP, BJ);     % Block Jacobi smoothing
    %xP = blkgs(xP, BJ, @(xP) bP - AP*xP); % Block GS smoothing
end
% -------------
% Visualization
% -------------
%x(perm) = xP;
%x(perm) = Pt(xP); % deflation
%x(perm) = V'*(VAV\(V*bP)) + Pt(xP); % deflation
%mesh(reshape(x, nx, nx))

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
