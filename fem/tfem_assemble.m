function [A,b]=tfem_assemble(coo,con,f)
n=size(coo,1); % number of vertices
ne=size(con,1); % number of elements (triangles)
c=3e8;
k=(2*pi*f)/c;
% A=spalloc(n,n,7*n); % approximate number of nonzeros is used here
indi=zeros(7*n,1);
indj=zeros(7*n,1);
vval=zeros(7*n,1);
cnt=1;
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
            indi(cnt) = con(e,i);
            indj(cnt) = con(e,j);
            vval(cnt) = -W(i,j)+(k^2)*M(i,j);
            cnt = cnt + 1;
%             A(con(e,i),con(e,j))=A(con(e,i),con(e,j))-W(i,j)+(k^2)*M(i,j);
        end
        b(con(e,i))=b(con(e,i))+f(i);
    end
end
A = sparse(indi,indj,vval,n,n);
end