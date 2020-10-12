function M = spai(A, lfill)
%SPAI Sparse Approximate Inverse.
%
%   Reference Paper:
%
%   Grote, M.J. and Huckle, T., 1997.
%   Parallel preconditioning with sparse approximate inverses.
%   SIAM Journal on Scientific Computing, 18(3), pp.838-853.
%
%   https://epubs.siam.org/doi/abs/10.1137/S1064827594276552?journalCode=sjoce3
%

n = length(A);
I = speye(n);
G = A^lfill; % inefficient, use BFS instead
M = G;
for i = 1:n
    ind = find(G(i,:)); % TODO: check if rows or cols
%     k = [];
%     for j = 1:length(ind)
%         k = [k find(A(:,ind(j)))'];
%     end
    [k, ~] = find(A(:,ind));
    k = unique(k);
    AA = A(k,ind);
    M(ind,i) = (AA'*AA) \ (AA'*I(k,i));
%     M(ind,i) = AA \ I(k,i);
%     spy(A(k,ind))
end
end
