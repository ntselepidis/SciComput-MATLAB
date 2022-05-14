function x = conjgrad(A, b, tol, maxit, M, dummy, x0)
    flag = 0;
    normb = norm(b);
    tolb = tol*normb;
    if nargin < 7
        x = zeros(length(b), 1);
    else
        x = x0;
    end
    r = b - A(x);
    z = M(r);
    p = z;
    gamma = r'*z;
    for i = 1 : maxit
        Ap = A(p);
        alpha = gamma / (p'*Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;
        normr = norm(r);
        if normr < tolb
            flag = 1;
            disp(['conjgrad converged in ' num2str(i) ' iterations with relative residual ' num2str(normr/normb) '.'])
            break;
        end
        z = M(r);
        gamma_old = gamma;
        gamma = r'*z;
        beta = gamma / gamma_old;
        p = z + beta*p;
    end
    if flag == 0
        disp(['conjgrad failed to converge within ' num2str(maxit) ' iterations.'])
    end
end
