clear;clc;close all;%rng(1);
a=0;
b=1;

f=0*5e8;    % if f=0 then Helmholtz -> Poisson

m=64;       % number of subintervals per dimension
ndoms=8;    % number of subdomains
fem_type=0; % 0 -> triangular elements, 1 -> square elements
lshape=1;   % 0 -> square domain, 1-> L-shape domain

% discretize domain
if ( fem_type )
    [h,ne,n,coo,con,bounds]=qfem(a,b,m);
else
    [h,ne,n,coo,con,bounds]=tfem(a,b,m);
end

if ( lshape )
    [coo,con,bounds,~,~]=lshape_fem(a,b,m,coo,con);
    ne=size(con,1);
end

% compute element centers and partition domain geometrically
ce=zeros(ne,2);
for i=1:ne
    xc=sum( coo( con(i,:), 1 ) ) / size(con,2);
    yc=sum( coo( con(i,:), 2 ) ) / size(con,2);
    ce(i,:)=[xc yc];
end
part=kmeans(ce,ndoms); % domain partitioning

% plot partitioning
clr=['b' 'r' 'g' 'c' 'm' 'y']';
clr=repmat(clr,ceil(ndoms/length(clr)),1);
figure, hold on;
for i=1:ndoms
    ind=(part==i);
    plot(ce(ind,1),ce(ind,2),strcat(clr(i),'*'));
end
hold off;

elems=(1:ne)';        % element indices
ielem=cell(ndoms,1);  % element indices of each domain
nelem=zeros(ndoms,1); % element number per domain
dcon1=cell(ndoms,1);  % domain connectivity matrix    (global numbering)
dcon2=cell(ndoms,1);  % domain connectivity matrix    ( local numbering)
All=cell(ndoms,1);    % vertex indices of each domain (global numbering) 
Ac=cell(ndoms,1);     % subdomain stiffness matrices 
bc=cell(ndoms,1);     % subdomain rhs vectors
xc=cell(ndoms,1);     % subdomain solution vectors
blk=zeros(ndoms+1,1); % subdomain index bounds
blk(1)=1;

% assemble subdomain stiffness matrices
for i=1:ndoms
    ielem{i}=elems(part==i);
    nelem(i)=length(ielem{i});
    dcon1{i}=con(ielem{i},:);
    
    All{i}=unique( dcon1{i} );
    mapp=zeros( length(coo), 1 );
    mapp(All{i})=1:length(All{i});
    dcon2{i}=mapp( dcon1{i} );
    
    nv=length(All{i});
    dcoo=zeros(nv,2);
    for j=1:nv
        dcoo(j,:)=coo(All{i}(j),:);
    end
    
    if ( fem_type )
        [Ac{i},bc{i}]=qfem_assemble(dcoo,dcon2{i},f);
    else
        [Ac{i},bc{i}]=tfem_assemble(dcoo,dcon2{i},f);
    end
    
    % remove boundaries ( dirichlet (0) )
    dbnds=mapp(bounds);
    dbnds=dbnds(dbnds>0);
    Ac{i}(dbnds,:)=[];
    Ac{i}(:,dbnds)=[];
    bc{i}(dbnds)=[];
    All{i}=setdiff(All{i},bounds);
    blk(i+1)=blk(i)+length(All{i});
end

p=vertcat(All{:});
n=size(coo,1);
deg=zeros(n,1);
for i=1:ndoms
    deg(All{i})=deg(All{i})+1;
end
s=find(deg>1);
deg=deg(deg>1);

% serial B and Bd computation
B=spalloc(sum(deg-1),length(p),sum( 2*(deg-1) ));
idx=1;
for i=1:length(s)
    ind=find(p==s(i));
    for j=1:deg(i)-1
        B(idx,ind(j))=1;
        B(idx,ind(j+1))=-1;
        idx=idx+1;
    end
end

Bc = cell(ndoms,1);
for i=1:ndoms
    Bc{i} = B(:,blk(i):blk(i+1)-1);
end

nl=size(B,1);
S=zeros(nl,nl);
g=zeros(nl,1);
for i=1:ndoms
    S = S + Bc{i} * ( Ac{i} \ Bc{i}' );
    g = g + Bc{i} * ( Ac{i} \ bc{i}  );
end
% AA=blkdiag(Ac{:});
% lambda=pcg(-S,-g,1e-10,100);
% lambda=pcg(-S,-g,1e-10,100,@(y) -B*(AA*(B'*y)));
lambda = S\g;
% lambda = 0*lambda;
for i=1:ndoms
    xc{i} = Ac{i} \ ( bc{i} - Bc{i}'*lambda ); 
end

x=zeros(n,1);
for i=1:ndoms
    x(All{i}) = xc{i};
end

figure, trimesh(con,coo(:,1),coo(:,2),x), % plot solution
title(sprintf('Solution of Helmholtz PDE for f = %e',f));