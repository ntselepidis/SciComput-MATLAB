clear;clc;close all;
a=-1;
b=1;
m=64;            % number of subintervals per dimension
h=(b-a)/m;       % mesh size
ne=m^2;          % number of elements (squares)
n=(m+1)^2;       % number of vertices
coo=zeros(n,2);  % coordinate matrix
con=zeros(ne,4); % connectivity matrix
t=(a:h:b)';
for i=1:(m+1)
    idx=(i-1)*(m+1)+1:i*(m+1);
    coo(idx,1)=t;
    coo(idx,2)=a+(i-1)*h;
end
idx=1;
for i=1:m
    for j=1:m
        con(idx,:)=[j j+1 j+m+2 j+m+1]+(i-1)*(m+1);
        idx=idx+1;
    end
end
south=1:m+1;
west=(m+2):(m+1):((m+1)^2-2*m-1);
north=south+(m+1)*m;
east=west+m;
% bounds=[south north west east]; % Dirichlet (0) on all bounds
bounds=[south (m+2):(m+1):((m+1)^2-m)]; % Dirichlet (0) on SW bounds
                                        % Neumann (0) elsewhere
f=5e8;
c=3e8;
k=(2*pi*f)/c;
A=spalloc(n,n,(3*m-2)^2);
b=zeros(n,1);
W=zeros(4);
M=zeros(4);
f=zeros(4,1);
for e=1:ne
    x1=coo(con(e,1),1);
    x2=coo(con(e,2),1);
    x3=coo(con(e,3),1);
    x4=coo(con(e,4),1);
    y1=coo(con(e,1),2);
    y2=coo(con(e,2),2);
    y3=coo(con(e,3),2);
    y4=coo(con(e,4),2);
    y41=y4-y1;
    y12=y1-y2;
    x14=x1-x4;
    x21=x2-x1;
    J=0.25*((x2-x1)*(y4-y1)-(x4-x1)*(y2-y1));
    W(1,1)=(2*x14^2+3*x14*x21+2*x21^2+2*y12^2+3*y12*y41+2*y41^2)/(24*J);
    W(2,2)=(2*x14^2-3*x14*x21+2*x21^2+2*y12^2-3*y12*y41+2*y41^2)/(24*J);
    W(3,3)=(2*x14^2+3*x14*x21+2*x21^2+2*y12^2+3*y12*y41+2*y41^2)/(24*J);
    W(4,4)=(2*x14^2-3*x14*x21+2*x21^2+2*y12^2-3*y12*y41+2*y41^2)/(24*J);
    W(1,2)=-(2*x14^2-x21^2-y12^2+2*y41^2)/(24*J);
    W(1,3)=-(x14^2+3*x14*x21+x21^2+y12^2+3*y12*y41+y41^2)/(24*J);
    W(1,4)=(x14^2-2*x21^2-2*y12^2+y41^2)/(24*J);
    W(2,3)=(x14^2-2*x21^2-2*y12^2+y41^2)/(24*J);
    W(2,4)=-(x14^2-3*x14*x21+x21^2+y12^2-3*y12*y41+y41^2)/(24*J);
    W(3,4)=-(2*x14^2-x21^2-y12^2+2*y41^2)/(24*J);
    W(2,1)=W(1,2); W(3,1)=W(1,3); W(4,1)=W(1,4);
    W(3,2)=W(2,3); W(4,2)=W(2,4);
    W(4,3)=W(3,4);
    M(1,1)=(4*J)/9;
    M(2,2)=(4*J)/9;
    M(3,3)=(4*J)/9;
    M(4,4)=(4*J)/9;
    M(1,2)=(2*J)/9;
    M(1,3)=J/9;
    M(1,4)=(2*J)/9;
    M(2,3)=(2*J)/9;
    M(2,4)=J/9;
    M(3,4)=(2*J)/9;
    M(2,1)=M(1,2); M(3,1)=M(1,3); M(4,1)=M(1,4);
    M(3,2)=M(2,3); M(4,2)=M(2,4);
    M(4,3)=M(3,4);
    f(1)=J;
    f(2)=J;
    f(3)=J;
    f(4)=J;
    for i=1:4
        for j=1:4
            A(con(e,i),con(e,j))=A(con(e,i),con(e,j))-W(i,j)+(k^2)*M(i,j);
        end
        b(con(e,i))=b(con(e,i))+f(i);
    end
end
A(bounds,:)=[];
A(:,bounds)=[];
b(bounds)=[];
x=A\b;
in=setdiff(1:n,bounds);
out=bounds;
xx=zeros(n,1);
xx(in)=x;
xx(out)=0;
figure, trimesh(con,coo(:,1),coo(:,2),xx);
