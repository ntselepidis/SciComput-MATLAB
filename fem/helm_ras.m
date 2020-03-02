% addpath('../meshpart');warning off;
clear;clc;close all;rng(0);
a=-1;
b=+1;
m=256;  % number of subintervals per dimension
f=0*5e8; % if f=0 then Helmholtz -> Poisson

ndoms  = 256;
nbasis = 4;

fem_type=0; % choose 0 for triangles or 1 for squares

if (fem_type)
    [h,ne,n,coo,con,bounds]=qfem(a,b,m);
    [A,b]=qfem_assemble(coo,con,f);
else
    [h,ne,n,coo,con,bounds]=tfem(a,b,m);
    [A,b]=tfem_assemble(coo,con,f);
end

bounds=[1:(m+1) (m+2):(m+1):((m+1)^2-m)]; % Dirichlet(0) on NE

in=setdiff(1:n,bounds);
xx=zeros(n,1);
A(bounds,:)=[];
A(:,bounds)=[];
b(bounds)=[];

% Graph Partitioning
G = 0.5*(abs(A)+abs(A)'); % undirected graph
map = PartSparseMat(G,ndoms);
% map = specdice(G,log2(ndoms))+1; 
ndoms = max(map);

% Classify Inner-Outer
[In, Out, All, blk] = inner_outer(G, map);

% RAS Preprocessing 
RAS.All = All;
[RAS.All2, ~] = form_overlap(G, map, All, Out);
[RAS.L, RAS.U, RAS.P, RAS.Q] = factorize(A, RAS.All2);

% SHEM Coarse Space for RAS
[R, basis] = shem(A, In, Out, All, nbasis);
RAR = R*A*R';

% Solve using RAS-BiCGSTAB
ex = bicgstab( A, b, 1e-8, 500, @(y) rasprec(y,RAS) );
ex = bicgstab( A, b, 1e-8, 500, @(y) rasprec(y,RAS) + R'*(RAR\(R*y)) );

xx(in)=ex;
xx(bounds)=0;

figure, trimesh(con,coo(:,1),coo(:,2),xx), % plot solution
title(sprintf('Solution of Helmholtz PDE for f = %e',f));