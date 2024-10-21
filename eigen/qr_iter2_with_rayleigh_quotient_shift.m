function [Qp,T]=qr_iter2_with_rayleigh_quotient_shift(H,tol,maxit)
[n,~]=size(H);
Qp=eye(n,n); % Qp = Q(1)*Q(2)*...*Q(maxit).
k=0;
iters=-1;
for m=n:-1:2
    while (1)
        k = k + 1;
        sigma_k = H(m,m);
        [Q,H] = hess_qr(H - sigma_k*eye(n,n)); % H = QR, Hbar = RQ, H = Hbar, in-place
        H = H + sigma_k*eye(n,n);
        Qp = Qp*Q;
        if ( norm(H(m,m-1))<tol )
            iters=k;
            break;
        end
    end
end
if (iters~=-1)
    disp(['QR algorithm converged in ' num2str(iters) ' iterations.'])
else
    disp(['QR algorithm did not converge within ' num2str(maxit) ...
        ' iterations.'])
end
T = diag(H);
end
