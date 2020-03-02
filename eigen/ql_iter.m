function [Qp,D_new]=ql_iter(A,tol,maxit)
[n,k]=size(A);
perm=k:-1:1; % reverse column ordering -> compute QL instead of QR decomp.
Q=zeros(n,k);
Qp=eye(n,k); % Qp = Q(1)*Q(2)*...*Q(maxit).
D_new=diag(A);
iters=-1;
for i=1:maxit
    Q(:,perm)=gs(A(:,perm)); % Gram-Schmidt on columns of A(i) (reversely).
    L=Q'*A; % Compute lower triangular matrix L.
    Qp=Qp*Q;
    A=L*Q; % Similarity transformation: A(i+1)=L(i)*Q(i)=Q(i)'*A(i)*Q(i).
    D_old=D_new;
    D_new=diag(A);
    if ( norm(D_new-D_old)<tol )
        iters=i;
        break;
    end
end
if (iters~=-1)
    disp(['QL algorithm converged in ' num2str(iters) ' iterations.'])
else
    disp(['QL algorithm did not converge within ' num2str(maxit) ...
        ' iterations.'])
end
end