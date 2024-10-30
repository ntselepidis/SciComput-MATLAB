function [U,H]=hess_qr(H)
[m,n]=size(H);
c = zeros(n-1, 1);
s = zeros(n-1, 1);
G = cell(n-1,1);
for k=1:n-1
    [c(k), s(k)] = givens(H(k,k), H(k+1,k));
    G{k} = [c(k) s(k); -s(k)' c(k)];
    H(k:k+1,k:n) = G{k} * H(k:k+1,k:n);
end
for k=1:n-1
    H(1:k+1,k:k+1) = H(1:k+1,k:k+1) * G{k}';
end
U = eye(n,n);
for k=1:n-1
    U(1:n,k:k+1) = U(1:n,k:k+1) * [c(k), -s(k)'; s(k), c(k)];
end
end
