clear; clc; close all; rng(0);

ndoms = 8;
lambda = 7;

a  = 0;
b  = 1;
nx = 100;
n  = nx^2;
h  = (b-a)/(nx+1); 

L = -gallery('poisson',nx)/(h^2);

[y,x] = meshgrid( (a+h):h:(b-h) );
x = x(:);
y = y(:);

dx   = x - x.^2;
dy   = y - y.^2;
dxdy = dx.*dy;

f = 2*(dx+dy) + lambda*dxdy.*exp(dxdy);

F  = @(u) L*u + lambda*exp(u) - f;
DF = @(u) L   + spdiags( lambda*exp(u), 0, n, n );

% Graph Partitioning
G = 0.5*(abs(L) + abs(L)');
map = PartSparseMat(G,ndoms);
ndoms = max(map);

% MTI Setup based on the sparsity pattern of the Jacobian matrix

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

% Krylov Method Setup
KSM.type = 0; % Bi-CGSTAB
KSM.tol = 1e-8;
KSM.maxit = 100;
KSM.pc = MTI;
KSM.pcfun = @(x,pc) mti_prec(x,pc);

% Newton Method Setup
uinit = zeros(n,1);
tol   = 1e-8;
maxit = 20;

% u = newton(F, DF, uinit, tol, maxit);
u = newton_krylov(F, DF, uinit, tol, maxit, KSM);

mesh( reshape(x,nx,nx), reshape(y,nx,nx), reshape(u,nx,nx) );




