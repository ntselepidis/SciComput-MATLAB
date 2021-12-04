clear; clc; close all; rng(0);

% Number of grid-points per dimension, mesh size, and number of unknowns
nx = 200;
h  = 1/(nx+1);
n  = nx^2;

% Setup linear system
A = gallery('poisson',nx)/(h^2);
ex = (1:n)'/n;
b = A*ex;

% Number of subdomains
ndoms = 8;

% Graph Partitioning
G = 0.5*(abs(A) + abs(A)');
map = PartSparseMat(G,ndoms);
ndoms = max(map);

% MTI Setup
[In, Out, All, blk] = inner_outer(G, map);

[nbr_doms, T, couples] = neighbor_domains(G, map, Out);

sep = vertex_separators(G, map, Out, couples);

[In, All, blk] = update_dpnts(All, sep, length(G));

[s, deg, p] = global_separators(All, sep);

[B, Bd] = assemble_B(s, deg, p);

Bc  = cell(ndoms,1);
Bdc = cell(ndoms,1);
for i=1:ndoms
    idx = blk(i):blk(i+1)-1;
    Bc{i}  = B(:,idx);
    Bdc{i} = Bd(:,idx);
end

[MTI.All, MTI.Bc, MTI.Bdc] = deal(All, Bc, Bdc);

[MTI.L, MTI.U, MTI.P, MTI.Q] = factorize( A, MTI.All );

% Krylov Method Setup
KSM.type = 0; % Bi-CGSTAB
KSM.tol = 1e-8;
KSM.maxit = 500;

% Unpreconditoned KSM
x = ksm(A, b, KSM.tol, KSM.maxit, @(y) y, KSM.type);

% Preconditioned KSM with MTI
x = ksm(A, b, KSM.tol, KSM.maxit, @(y) mti_prec(y, MTI), KSM.type);

% Visualize solution
% mesh( reshape(x,nx,nx) );
