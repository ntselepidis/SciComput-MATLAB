clear; clc; close all; rng(0);

ndoms  = 256;
nbasis = 4;
gptype = 0; % 0 for kmeans | 1 for metis

tol    = 1e-8;
maxit  = 500;

x0 = 0;
x1 = 0.25;
m = 128;       % number of subintervals per dimension

E  = 210*1e6; % Elastic modulus (Young's modulus)
NU = 0.3;     % Poisson's ratio
t  = 0.025;   % thickness
w  = 100;     % distributed load
P  = 12.5;    % concentrated load
p  = 1;       % plane stress

[h,ne,n,coo,con,bounds] = tfem(x0,x1,m);
[coo,con,bounds,sep,dpnts] = lshape_fem_2(x0,x1,m,coo,con);
ne = size(con,1);
n  = size(coo,1);

% Vertex numbering per square subdomain
n1 = dpnts{1};
n2 = dpnts{2};
n3 = dpnts{3};

% Dirichlet Boundaries (0)
west1 = n1(1:(m+1):(m+1)^2-m);
west3 = n3(1+(m+1):(m+1):(m+1)^2-m);
bounds = [west1 west3];
bounds = [2*bounds-1 2*bounds]; % both ux and uy components are zero

% Form forces (rhs vector)
b = zeros(2*n,1);
north2 = n2( (1:(m+1)) + (m+1)*m );
b(2*north2) = -(w*x1)/(m+1);      % distributed load
b(2*n2(end)) = b(2*n2(end)) + P;  % concentrated load

% A = spalloc(2*n,2*n,12*n);
% Coordinate format of global stiffness matrix
indi = zeros(12*n,1);
indj = zeros(12*n,1);
vval = zeros(12*n,1);
cnt  = 1;

% assemble global stiffness matrix 
for e=1:ne
    x1 = coo(con(e,1),1);
    y1 = coo(con(e,1),2);
    x2 = coo(con(e,2),1);
    y2 = coo(con(e,2),2);
    x4 = coo(con(e,3),1);
    y4 = coo(con(e,3),2);
    M = LinearTriangleElementStiffness(E,NU,t,x1,y1,x2,y2,x4,y4,p);    
    for i=1:3
        for j=1:3
            ii=con(e,i);
            jj=con(e,j);
            indi(cnt) = 2*ii-1;
            indj(cnt) = 2*jj-1;
            vval(cnt) = M(2*i-1,2*j-1);
            cnt = cnt + 1;
            indi(cnt) = 2*ii-1;
            indj(cnt) = 2*jj;
            vval(cnt) = M(2*i-1,2*j);
            cnt = cnt + 1;
            indi(cnt) = 2*ii;
            indj(cnt) = 2*jj-1;
            vval(cnt) = M(2*i,2*j-1);
            cnt = cnt + 1;
            indi(cnt) = 2*ii;
            indj(cnt) = 2*jj;
            vval(cnt) = M(2*i,2*j);
            cnt = cnt + 1;
%             A(2*ii-1,2*jj-1) = A(2*ii-1,2*jj-1) + M(2*i-1,2*j-1);
%             A(2*ii-1,2*jj  ) = A(2*ii-1,2*jj  ) + M(2*i-1,2*j  );
%             A(2*ii  ,2*jj-1) = A(2*ii  ,2*jj-1) + M(2*i  ,2*j-1);
%             A(2*ii  ,2*jj  ) = A(2*ii  ,2*jj  ) + M(2*i  ,2*j  ); 
        end
    end
end

% A = sparse(indi,indj,vval,2*n,2*n);
A = sparse(indi(1:cnt-1),indj(1:cnt-1),vval(1:cnt-1),2*n,2*n);

disp('Matrix Assembly - Done');

ind = setdiff(1:length(A),bounds); % all dofs appart from boundaries

bnd = [west1 west3];    % bound points
inn = setdiff(1:n,bnd); % inner points

u = zeros(length(A),1); % allocate solution vector and fill with zeros

% ------------------------------------------------------------------------
% ----------------------- Solve using RAS-BiCGSTAB -----------------------
% ------------------------------------------------------------------------
AA = A(ind,ind);
bb = b(ind);

G = 0.5*( abs(AA) + abs(AA)' );

if ( gptype )
    % Metis Partitioning
    map = PartSparseMat( G, ndoms ); 
else
    % Geometric Partitioning
    part = kmeans(coo(inn,:),ndoms,'MaxIter',200);
    map = zeros(1,length(ind));
    for i=1:ndoms
        idx = find( part==i );
        map(2*idx-1) = i;
        map(2*idx)   = i;
    end
end
ndoms = max(map);

disp('Partitioning - Done');

[In, Out, RAS.All, blk]      = inner_outer( G, map );
[RAS.All2, ~]                = form_overlap( G, map, RAS.All, Out );
[RAS.L, RAS.U, RAS.P, RAS.Q] = factorize( AA, RAS.All2 );

% SHEM Coarse Space for RAS
[R, basis] = shem( AA, In, Out, RAS.All, nbasis );
RAR = R*AA*R';
[cL,cU,cP,cQ] = lu(RAR, 'vector'); 

disp('RAS Preprocessing - Done');

u(ind) = bicgstab( AA, bb, tol, maxit, @(y) rasprec(y,RAS) );
u(ind) = bicgstab( AA, bb, tol, maxit, ...
    @(y) rasprec(y,RAS) + R'*sp_solve(cL,cU,cP,cQ,R*y) );

% ------------------------------------------------------------------------

ux = u( ind(1:2:end) ); % x component of u (active dofs)
uy = u( ind(2:2:end) ); % y component of u (active dofs)

% bnd = [west1 west3];    % bound points
% inn = setdiff(1:n,bnd); % inner points

f = A*u;         % compute forces
fx = f(2*bnd-1); % x component of f on bounds
fy = f(2*bnd);   % y component of f on bounds

% plot solution
inc = 0.05;
figure, hold on, 
quiver(coo(inn,1), coo(inn,2), ux, uy, 'b');
quiver(coo(bnd,1), coo(bnd,2), fx, fy, 'r');
axis([x0-inc 2*x1+inc x0-inc 2*x1+inc]);
axis square; grid; xlabel('x'); ylabel('y');
title('displacements and bound reactions');
legend('displacements','bound reactions');
