function [A, b] = qfem_assemble(coo, con, f)
n = size(coo,1);  % number of vertices
ne = size(con,1); % number of elements (squares)
c = 3e8;
k = (2*pi*f)/c;
% A = spalloc(n,n,9*n); % approximate number of nonzeros is used here
indi = zeros(9*n,1);
indj = zeros(9*n,1);
vval = zeros(9*n,1);
cnt = 1;
b = zeros(n,1);
for e = 1:ne
    x1 = coo(con(e,1),1);
    x2 = coo(con(e,2),1);
    x3 = coo(con(e,3),1);
    x4 = coo(con(e,4),1);
    y1 = coo(con(e,1),2);
    y2 = coo(con(e,2),2);
    y3 = coo(con(e,3),2);
    y4 = coo(con(e,4),2);
    [W, M, f] = qfem_stiffness(x1, y1, x2, y2, x3, y3, x4, y4);
    for i = 1:4
        for j = 1:4
            indi(cnt) = con(e,i);
            indj(cnt) = con(e,j);
            vval(cnt) = -W(i,j) + (k^2)*M(i,j);
            cnt = cnt + 1;
%             A(con(e,i),con(e,j)) = A(con(e,i),con(e,j)) - W(i,j) + (k^2)*M(i,j);
        end
        b(con(e,i)) = b(con(e,i)) + f(i);
    end
end
A = sparse(indi, indj, vval, n, n);
end
