function [h,ne,n,coo,con,bounds]=qfem(a,b,m)
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
bounds=[south north west east];
end