function x = richardson(A, b, tol, maxit, M, dummy, x0)
    flag = 0;
    % initial estimate
    if nargin < 7
        x = zeros(length(b), 1);
    else
        x = x0;
    end
    % richardson iteration
    normb = norm(b);
    for i = 1 : maxit
        r = b - A(x);
        relres = norm(r) / normb;
        if (relres <= tol)
            disp(['richardson iteration converged in ' num2str(i) ' iterations with relative residual ' num2str(relres) '.'])
            flag = 1;
            break;
        end
        rhat = M(r);
        %omega = (rhat'*rhat) / (rhat'*(A(M(rhat))));
        omega = (r'*rhat) / (rhat'*A(rhat));
        x = x + omega * rhat;
    end
    if flag == 0
        disp(['richardson iteration failed to converge within ' num2str(maxit) ' iterations.'])
    end
end
