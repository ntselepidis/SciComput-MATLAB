function [A,b]=qfem_assemble(coo,con,f)
n=size(coo,1);  % number of vertices
ne=size(con,1); % number of elements (squares)
c=3e8;
k=(2*pi*f)/c;
% A=spalloc(n,n,9*n); % approximate number of nonzeros is used here
indi=zeros(9*n,1);
indj=zeros(9*n,1);
vval=zeros(9*n,1);
cnt=1;
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