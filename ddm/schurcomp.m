clear; clc; close all; rng(0);

nx    = 200;
ndoms = 32;
KSM   = 0;
tol   = 1e-8;
maxit = 500*(KSM==0) + 100*(KSM==1) + 500*(KSM==2);
delta = 1; % ATTENTION : CURRENT SCHURAS IMPLEMENTATION SUPPORTS ONLY ONE-LEVEL OVERLAP

A = gallery('poisson',nx);
n = length(A);
ex = (1:n)'/n;
b = A*ex;

% Graph Partitioning
G = 0.5*(abs(A)+abs(A)'); % undirected graph
map = PartSparseMat(G, ndoms);
ndoms = max(map);

% Classify subdomain components into inner and outer points 
[In, Out, All, blk] = inner_outer(G, map);

% --------------------------------------------------------------------- %
% ---------------------- Schur Complement Method ---------------------- %
% --------------------------------------------------------------------- %

% Schur Complement Method Preprocessing
SC = sc_prep(A, In, Out);

% Schur Complement Block Jacobi Setup
SBJ = scbj_prep(A, SC);

[out, L, U, P, Q, Aio, Aoo] = deal(SC.out, SC.L, SC.U, SC.P, SC.Q, SC.Aio, SC.Aoo);

% --------------------------------------------------------------------- %
% ---------------------- Overlap Visualization ------------------------ %
% --------------------------------------------------------------------- %

% delta = input('delta = ');
% 
% overlap = form_overlap_ml(G,In,Out,delta);
% 
% for i=1:ndoms
%     ovlp = setdiff( overlap{i}, All{i} );
%     RAS.All2{i} = [All{i} ovlp];
% end
% 
% nx = sqrt(n);
% h = 1 / (nx+1);
% [y, x] = meshgrid( h : h : 1-h );
% y = y(:); x = x(:);
% clr = ['b' 'r' 'k' 'm'];
% smb = ['s' 's' 's' 's'];
% figure, hold on;
% for i=1:ndoms
%     plot( x(All{i}), y(All{i}), strcat(clr(i),'*') );
% end
% for i=1:ndoms
%     plot( x(overlap{i}), y(overlap{i}), strcat(clr(5-i),smb(i)) );
% end
% hold off;
% 
% return

assert( delta == 1 );

overlap = form_overlap_ml( G, In, Out, delta );

ilu_setup.type = 'ilutp';
ilu_setup.droptol = 1e-3;
ilu_setup.milu = 'off';
ilu_setup.udiag = 0;
ilu_setup.thresh = 1;

s_droptol = 1e-2;

AS = scas_prep(A, SC, overlap, s_droptol, ilu_setup);

x = zeros(n,1);

% Preconditioned KSM Solvers for Interface Problem (Schur Complement)
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) speye(length(y))*y, KSM); % No Prec
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) blkjac(y,SBJ), KSM); % SBJ
x(out) = ksm(@(x) Smultx(x,SC), Srhs(b,SC), tol, maxit, @(y) schuras(y,AS), KSM); % SBJ
return
figure, mesh(reshape(x,nx,nx))

for j=1:ndoms
    x(In{j}) = sp_solve( L{j},U{j},P{j},Q{j}, b(In{j})-Aio{j}*x(out) );
end

% figure, mesh(reshape(x,nx,nx))
