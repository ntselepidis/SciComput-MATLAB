function [x, val_record] = richardson(A, b, tol, maxit, M, dummy, x0)
    flag = 0;
    normb = norm(b);
    tolb = tol*normb;
    if nargin < 7
        x = zeros(length(b), 1);
    else
        x = x0;
    end
    val1 = 0.5*( x'*A(x) ) - x'*b;
    val_record = zeros(maxit+1, 1);
    val_record(1) = val1;
    r = b - A(x);
    p = M(r);
    for i = 1 : maxit
        Ap = A(p);
        %alpha = (p'*p) / (p'*(A(M(p))));
        %alpha = (p'*p) / (p'*(M(Ap)));
        alpha = (r'*p) / (p'*Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;
        val0 = val1;
        val1 = 0.5*( x'*A(x) ) - x'*b;
        val_record(i+1) = val1;
        normr = norm(r);
        %if (normr <= tolb)
        %%if (abs(val1 - val0)/abs(val0) <= tol)
        %    flag = 1;
        %    disp(['richardson iteration converged in ' num2str(i) ' iterations with relative residual ' num2str(normr/normb) '.'])
        %    break;
        %end
        p = M(r);
    end
    %if flag == 0
    %    disp(['richardson iteration failed to converge within ' num2str(maxit) ' iterations.'])
    %end
end
