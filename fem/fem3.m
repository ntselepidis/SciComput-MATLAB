clear;clc;close all;
a=-1;
b=1;
m=64;            % number of subintervals per dimension
h=(b-a)/m;       % mesh size
ne=2*m^2;        % number of elements (triangles)
n=(m+1)^2;       % number of vertices
coo=zeros(n,2);  % coordinate matrix
con=zeros(ne,3); % connectivity matrix
t=(a:h:b)';
for i=1:(m+1)
    idx=(i-1)*(m+1)+1:i*(m+1);
    coo(idx,1)=t;
    coo(idx,2)=a+(i-1)*h;
end
idx=1;
for i=1:m
    for j=1:m
        con(idx,:)=[j j+1 j+m+1]+(i-1)*(m+1);
        con(idx+1,:)=[j+m+2 j+m+1 j+1]+(i-1)*(m+1);
        idx=idx+2;
    end
end
south=1:m+1;
west=(m+2):(m+1):((m+1)^2-2*m-1);
north=south+(m+1)*m;
east=west+m;
% bounds=[south north west east];
bounds=[south (m+2):(m+1):(m+1)^2-m];

f=5e8;
c=3e8;
k=(2*pi*f)/c;
A=spalloc(n,n,m*(3*m-2)+(m-1)*(4*m-2));
b=zeros(n,1);
W=zeros(3);
M=zeros(3);
f=zeros(3,1);
for e=1:ne
    x1=coo(con(e,1),1);
    x2=coo(con(e,2),1);
    x4=coo(con(e,3),1);
    y1=coo(con(e,1),2);
    y2=coo(con(e,2),2);
    y4=coo(con(e,3),2);
    y41=y4-y1;
    y12=y1-y2;
    x14=x1-x4;
    x21=x2-x1;
    J=0.25*((x2-x1)*(y4-y1)-(x4-x1)*(y2-y1));
    W(1,1)=(x14^2 + 2*x14*x21 + x21^2 + y12^2 + 2*y12*y41 + y41^2)/(8*J);
    W(2,2)=(x14^2 + y41^2)/(8*J);
    W(3,3)=(x21^2 + y12^2)/(8*J);
    W(1,2)=-(x14^2 + x21*x14 + y41^2 + y12*y41)/(8*J);
    W(1,3)=-(x21^2 + x14*x21 + y12^2 + y41*y12)/(8*J);
    W(2,3)=(x14*x21 + y12*y41)/(8*J);
    W(2,1)=W(1,2); W(3,1)=W(1,3);
    W(3,2)=W(2,3);
    M(1,1)=J/3;
    M(2,2)=J/3;
    M(3,3)=J/3;
    M(1,2)=J/6;
    M(1,3)=J/6;
    M(2,3)=J/6;
    M(2,1)=M(1,2); M(3,1)=M(1,3);
    M(3,2)=M(2,3);
    f(1)=(2*J)/3;
    f(2)=(2*J)/3;
    f(3)=(2*J)/3;
    for i=1:3
        for j=1:3
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