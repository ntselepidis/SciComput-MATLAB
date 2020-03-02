function M = probe2( A, d, n )
% n = length(A);
k = min( [2*d+1,n] );
if k == 1
    M = spdiags( feval(A,ones(n,1)), 0, n, n );
else
    v = rem( [1:n]'*ones(1,k), k ) == ones(n,1)*[1:(k-1),0];
    av = feval(A,v);
    M = spalloc(n,n,n);
    for c=1:k
        for i=c:k:n
            M( max([i-d 1]):min([i+d n]), i ) = av( max([i-d 1]):min([i+d n]), c );
        end
    end
end
end