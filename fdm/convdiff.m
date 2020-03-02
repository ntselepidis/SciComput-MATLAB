function A = convdiff(dim,n,alpha,epsilon)

assert( n > 0 );
assert( dim == 2 || dim == 3 );
assert( length(alpha) == dim );

h = 1/(n+1);
e = ones(n,1);
I = speye(n);
T = spdiags(e*[-1 2 -1],-1:1,n,n) / (h^2); % second derivative
C = spdiags(e*[-1 1],[-1 1],n,n) / (2*h);  % first derivative

if (dim==2)
    TT = kron(I,T) + kron(T,I);
    CC = alpha(1)*kron(I,C) +...
         alpha(2)*kron(C,I);
else
    TT = kron(kron(I,I),T) + kron(kron(I,T),I) + kron(kron(T,I),I);
    CC = alpha(1)*kron(kron(I,I),C) +...
         alpha(2)*kron(kron(I,C),I) +...
         alpha(3)*kron(kron(C,I),I);
end

A = epsilon*TT + CC;

end