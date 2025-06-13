function R = coarsening(A)
% clear; clc;
% A = gallery('poisson', 4);
% U = 0 UNMATCHED
% F = 1 FINE
% C = 2 COARSE
G = abs(A) + abs(A)';
n = length(A);
v = zeros(n, 1);

for i = 1:n
    if ( v(i) == 0 )
        v(i) = 2;
        neib = setdiff(find(G(:,i)), i);
        v(neib) = 1;
    end
end    

fine = find( v == 1 );
nf = length(fine);
coarse = find( v == 2 );
nc = length(coarse);
map(coarse) = (1:nc);
nzmax = 7*nc;
R = spalloc(nc, n, nzmax);
for i = 1:nc
    R(i,coarse(i)) = 1;
end
for i = 1:nf
    neib = setdiff( find(G(:,fine(i))), fine(i) );
    cneib = neib( v(neib) == 2 );
    deg = length(cneib);
    R(map(cneib),fine(i)) = 1/deg;
end
end
