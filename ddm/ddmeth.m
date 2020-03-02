clear; clc; close all; rng(0);

ndoms  = 256; % number of subdomains
prob   = 7;  % probing parameter
nbasis = 5;  % number of basis functions for SHEM coarse grid

KSM   = 2;   % 0 for BiCGSTAB | 1 for GMRES(20) | 2 for full GMRES
tol   = 1e-8;
maxit = 500*(KSM==0) + 100*(KSM==1) + 500*(KSM==2);

% A = convdiff(2,100,[1.0 1.0],1e-4);    % 2D Convection-Diffusion
% A = convdiff(3,25,[1.0 0.0 0.0],1e-4);    % 3D Convection-Diffusion
% A = convdiff3D(0,1,25,1,1,1e-3); % 3D Convection-Diffusion (Haidar)
% A = convdiff2D(0,1,200,1,1,1e-3); % 2D Convection-Diffusion (Haidar)

A = convdiff(2,256,[0.0 0.0],1);       % 2D Poisson
% A = convdiff(2,256,[2.5 0.0],0.1);     % 2D Convection-Diffusion
% A = convdiff(2,256,[1.0 0.0],1e-4);    % 2D Convection-Diffusion
% A = convdiff(3, 32,[0.0 0.0 0.0],1);   % 3D Poisson
% A = convdiff(3, 32,[2.5 2.5 2.5],0.1); % 3D Convection-Diffusion
% A = convdiff3D(0,1,32,1,1,1e-3); % 3D Convection-Diffusion (Haidar)
n = length(A);
b = ones(n,1);
% ex = (1:n)'/n; b = A*ex;

unsymm = nnz(A-A')/nnz(A)

% Graph Partitioning
G = 0.5*(abs(A)+abs(A)'); % undirected graph
map = PartSparseMat(G,ndoms);
ndoms = max(map);

% Classify subdomain components into inner and outer points 
[In, Out, All, blk] = inner_outer(G, map);

% --------------------------------------------------------------------- %
% ---------------------- Schur Complement Method ---------------------- %
% --------------------------------------------------------------------- %

% Schur Complement Method Preprocessing
SC = sc_prep(A,In,Out);

% Schur Complement Block Jacobi Setup
SBJ = scbj_prep(A,SC);

% Schur Complement Block Probing Setup
SBP = scbp_prep(A,SC,prob);

[out,L,U,P,Q,Aio,Aoo] = deal(SC.out,SC.L,SC.U,SC.P,SC.Q,SC.Aio,SC.Aoo);

% --------------------------------------------------------------------- %
% -------------------------- One-Level DDMs --------------------------- %
% --------------------------------------------------------------------- %

% Block Jacobi Preprocessing
BJ.blk = blk;
[BJ.L, BJ.U, BJ.P, BJ.Q] = factorize(A, All);

% RAS Preprocessing
RAS.All = All;
[RAS.All2, ~] = form_overlap(G, map, All, Out);
[RAS.L, RAS.U, RAS.P, RAS.Q] = factorize(A, RAS.All2);

% permute matrix A
perm = horzcat( All{:} )'; 
AP = A(perm,perm); 
bP = b(perm);

% --------------------------------------------------------------------- %
% ----------------- Coarse Spaces for Two-Level DDMs ------------------ %
% --------------------------------------------------------------------- %

% Nicolaides Coarse Space for block Jacobi
V = aggregate(AP, blk);
VAV = V*AP*V';

% Nicolaides Coarse Space for RAS
W = spalloc(ndoms,n,n);
for i=1:ndoms
    W(i,All{i}) = 1;
end
WAW = W*A*W';

% SHEM Coarse Space for RAS
[R, basis] = shem(A, In, Out, All, nbasis);
RAR = R*A*R';

% Nicolaides Coarse Space for Schur Complement Method
Z = spalloc( ndoms, length(out), length(out) );
mask = zeros(length(A),1);
mask(out) = 1:length(out);
for i=1:ndoms
    Z(i,mask(Out{i})) = 1;
end
ZAZ = Z*Aoo*Z';
ZSZ = ZAZ - Z*blkdiag( SBJ.T{:} )*Z';

% --------------------------------------------------------------------- %
% ------------------ Preconditioned Krylov Solvers -------------------- %
% --------------------------------------------------------------------- %

x = zeros(n,1);

% Preconditioned KSM Solvers for Interface Problem (Schur Complement)
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBJ), KSM); % SBJ
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBP), KSM); % SBP

x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBJ)+Z'*(ZAZ\(Z*y)), KSM); % SBJ + approx Nico
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBP)+Z'*(ZAZ\(Z*y)), KSM); % SBP + approx Nico

x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBJ)+Z'*(ZSZ\(Z*y)), KSM); % SBJ + exact Nico
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBP)+Z'*(ZSZ\(Z*y)), KSM); % SBP + exact Nico

for j=1:ndoms
    x(In{j}) = sp_solve( L{j},U{j},P{j},Q{j}, b(In{j})-Aio{j}*x(out) );
end

% Preconditioned KSM Solvers for Ax=b
x = ksm( AP, bP, tol, maxit, @(y) blkjac(y,BJ), KSM );                    % Block Jacobi
x = ksm( AP, bP, tol, maxit, @(y) blkjac(y,BJ) + V'*(VAV\(V*y)), KSM );   % Block Jacobi + Nico
x = ksm( A,  b,  tol, maxit, @(y) rasprec(y,RAS), KSM );                  % RAS
x = ksm( A,  b,  tol, maxit, @(y) rasprec(y,RAS) + W'*(WAW\(W*y)), KSM ); % RAS + Nico
x = ksm( A,  b,  tol, maxit, @(y) rasprec(y,RAS) + R'*(RAR\(R*y)), KSM ); % RAS + SHEM
