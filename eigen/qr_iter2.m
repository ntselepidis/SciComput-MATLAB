function [Qp,D_new]=qr_iter2(H,tol,maxit)
[n,k]=size(H);
Qp=eye(n,k); % Qp = Q(1)*Q(2)*...*Q(maxit).
D_new=diag(H);
iters=-1;
for i=1:maxit
    [Q,H] = hess_qr(H); % H = QR, Hbar = RQ, H = Hbar, in-place
    Qp = Qp*Q;
    D_old=D_new;
    D_new=diag(H);
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
