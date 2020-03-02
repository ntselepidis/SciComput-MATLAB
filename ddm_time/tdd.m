clear; clc; close all; rng(0);

alpha = 0.5;  % diffusion coefficient
epsil = 0;    % convection coefficient
nx    = 32;   % number of unknowns per dimension
nt    = 101;  % number of timesteps (total)
m     = 10;   % number of steps in each "grouped" timestep
ndoms = 8;    % number of subdomains
tol   = 1e-8; % KSM tolerance
maxit = 100;  % KSM max number of iterations
prob  = 50;   % probing parameter (set to -1 to turn probing OFF)

% space parameters
x0 = 0;
x1 = 1;
h  = (x1-x0)/(nx+1); % mesh size
n  = nx^2;

[y,x] = meshgrid( (x0+h):h:(x1-h) ); 
x = x(:);
y = y(:);

% time parameters
t0 = 0;
t1 = 0.10;
dt = (t1-t0)/(nt-1); % timestep size
nT = floor( nt / m ); % number of "grouped" timesteps

% single timestep matrix assembly
I = speye(n);
C = spdiags(ones(nx,1)*[-1 1],[-1 1],nx,nx);
A = I + ((alpha*dt)/(h^2))*gallery('poisson',nx) - ((epsil*dt)/(2*h))*kron(speye(nx),C); 

% multiple timestep matrix assembly
K = matrix_unroll(A, m);

% initial condition setup
f = zeros( length(K), 1 );
f(1:n) = sin(pi*x).*sin(pi*y);

% Graph Partitioning
G = 0.5*(abs(K)+abs(K)'); % undirected graph
map = PartSparseMat(G,ndoms);
ndoms = max(map);

% Classify subdomain components into inner and outer points 
[In, Out, All, blk] = inner_outer(G, map);

% Schur Complement Method Preprocessing
SC = sc_prep(K,In,Out);

tic
if (prob == -1)
    % Schur Complement Block Jacobi Setup
    SBJ = scbj_prep(K,SC);
    disp('Probing OFF');
else
    % Schur Complement Block Probing Setup
    SBP = scbp_prep(K,SC,prob);
    SBJ = SBP;
    disp('Probing ON');
end
toc

% Extract data from SC data structure
[In,out,L,U,P,Q,Kio] = deal(SC.In,SC.out,SC.L,SC.U,SC.P,SC.Q,SC.Aio);

% % Coarse Space Construction
% V = spalloc( ndoms, length(out), length(out) );
% mask = zeros(length(K),1);
% mask(out) = 1:length(out);
% for i=1:ndoms
%     V(i,mask(Out{i})) = 1;
% end
% VKV = V*SC.Aoo*V';
% VSV = VKV - V*blkdiag( SBJ.T{:} )*V';

% Start Simulation
umin = min(f);
umax = max(f);
x = reshape(x,nx,nx);
y = reshape(y,nx,nx);
mesh( x, y, reshape(f(1:n),nx,nx) ); % plot initial condition
axis( [x0+h x1-h x0+h x1-h umin umax] );

u = zeros( length(K), 1 );
for i=1:nT

    % Solve interface problem
    u(out) = pbicgstab2(@(x) Smultx(x,SC), Srhs(f,SC), tol, maxit, @(x) blkjac(x,SBJ));
%     u(out) = pbicgstab2(@(x) Smultx(x,SC), Srhs(f,SC), tol, maxit, @(x) blkjac(x,SBJ) + V'*(VKV\(V*x)) );
%     u(out) = pbicgstab2(@(x) Smultx(x,SC), Srhs(f,SC), tol, maxit, @(x) blkjac(x,SBJ) + V'*(VSV\(V*x)) );

    % Back substitute to inner dofs
    for j=1:ndoms
        u(In{j}) = sp_solve( L{j}, U{j}, P{j}, Q{j}, f(In{j}) - Kio{j} * u(out) );
    end
    
    f(1:n) = u(end-n+1:end);

    % Plot batch solution
    for j=1:m
        II = (j*n+1):(j+1)*n;
        mesh( x, y, reshape(u(II),nx,nx) );
        axis( [x0+h x1-h x0+h x1-h umin umax] );
        getframe;
    end
    
end
