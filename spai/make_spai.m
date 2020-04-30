function M = make_spai(A, lfill)
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
