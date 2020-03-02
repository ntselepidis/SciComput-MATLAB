function [Qp,D_new]=qr_iter(A,tol,maxit)
[n,k]=size(A);
Qp=eye(n,k); % Qp = Q(1)*Q(2)*...*Q(maxit).
D_new=diag(A);
iters=-1;
for i=1:maxit
    [Q,R]=house_qr(A);
    Qp=Qp*Q;
    A=R*Q; % Similarity transformation: A(i+1)=R(i)*Q(i)=Q(i)'*A(i)*Q(i).
    D_old=D_new;
    D_new=diag(A);
    if ( norm(D_new-D_old)<tol )
        iters=i;
        break;
    end
end
if (iters~=-1)
    disp(['QR algorithm converged in ' num2str(iters) ' iterations.'])
else
    disp(['QR algorithm did not converge within ' num2str(maxit) ...
        ' iterations.'])
end
end