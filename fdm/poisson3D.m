function K3D = poisson3D(n)
    K1D = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n); % 1d Poisson matrix
    I1D = speye(size(K1D));                      % 1d Identity matrix
    K2D = kron(K1D,I1D) + kron(I1D,K1D);         % 2d Poisson matrix
    I2D = speye(size(K2D));                      % 2d Identity matrix
    K3D = kron(K2D,I1D)+kron(I2D,K1D);           % 3d Poisson matrix
end