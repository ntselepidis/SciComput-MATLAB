function K = matrix_unroll(A,m)

n = length(A);
[Ai,Aj,Av] = find(A);

Kr  = (m+1)*n;
Knz = m*nnz(A)+(m+1)*n;
Ki  = zeros(Knz,1);
Kj  = zeros(Knz,1);
Kv  = zeros(Knz,1);

Ki(1:n) = 1:n;
Kj(1:n) = 1:n;
Kv(1:n) = 1;
% si = n;
ei = n;
for i=1:m
    si = ei + 1;
    ei = ei + n;
    Ki(si:ei) = ((1:n) + i*n);
    Kj(si:ei) = ((1:n) + i*n - n);
    Kv(si:ei) = -1;
end
for i=1:m
    si = ei + 1;
    ei = ei + nnz(A);
    Ki(si:ei) = (Ai + i*n);
    Kj(si:ei) = (Aj + i*n);
    Kv(si:ei) = Av;
end

K = sparse( Ki, Kj, Kv, Kr, Kr );

end