function [V, H] = arnoldi(A, v1, m)
  % ARNOLDI  m-step Arnoldi iteration for matrix A, start vector v1.
  % Returns:
  %   V : n x (m+1) matrix, Arnoldi basis vectors (orthonormal columns)
  %   H : (m+1) x m upper Hessenberg matrix

  n = length(v1);
  V = zeros(n, m+1);
  H = zeros(m+1, m);

  % Normalize starting vector
  V(:,1) = v1 / norm(v1);

  for j = 1:m
    % 1) Apply A
    w = A * V(:,j);

    % 2) Modified Gram–Schmidt against existing basis
    for i = 1:j
      H(i,j) = V(:,i)' * w;
      w = w - H(i,j) * V(:,i);
    end

    % 3) Next basis vector and subdiagonal entry
    H(j+1,j) = norm(w);
    if H(j+1,j) == 0
      % Breakdown: Krylov subspace has dimension j
      V = V(:,1:j);   % trim
      H = H(1:j+1,1:j);
      return;
    end
    V(:,j+1) = w / H(j+1,j);
  end
end

