function [Qp,T]=qr_iter2_with_rayleigh_quotient_shift(H,tol,maxit)
[n,~]=size(H);
Qp=eye(n,n); % Qp = Q(1)*Q(2)*...*Q(maxit).
k=0;
iters=0;
maxit_inner = floor(maxit / n);
%maxit_inner = n;
for m=n:-1:2
    for iter_inner = 1 : maxit_inner
        k = k + 1;
        sigma_k = H(m,m);
        [~,H(1:m,1:m)] = hess_qr(H(1:m,1:m) - sigma_k*eye(m,m)); % H = QR, Hbar = RQ, H = Hbar, in-place
        H(1:m,1:m) = H(1:m,1:m) + sigma_k*eye(m,m);
        %Qp = Qp*Q;
        if ( abs(H(m,m-1)/H(m,m))<tol )
            iters=iters+k;
            break;
        end
    end
end
if (iters~=0)
    disp(['QR algorithm converged in ' num2str(iters) ' iterations.'])
else
    disp(['QR algorithm did not converge within ' num2str(maxit) ...
        ' iterations.'])
end
T = diag(H);
end
