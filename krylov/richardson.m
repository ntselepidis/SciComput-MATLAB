function [x, val_record] = richardson(A, b, tol, maxit, M, dummy, x0)
    flag = 0;
    % initial estimate
    if nargin < 7
        x = zeros(length(b), 1);
    else
        x = x0;
    end
    % richardson iteration
    normb = norm(b);
    val0 = 0.5*( x'*A(x) ) - x'*b;
    val1 = inf;
    val_record = zeros(maxit+1, 1);
    val_record(1) = val0;
    for i = 1 : maxit
        r = b - A(x);
        %relres = norm(r) / normb;
        %if (relres <= tol)
        %if (abs(val1 - val0)/abs(val0) <= tol)
        %    disp(['richardson iteration converged in ' num2str(i) ' iterations with relative residual ' num2str(relres) '.'])
        %    flag = 1;
        %    break;
        %end
        %val0 = val1;
        rhat = M(r);
        %omega = (rhat'*rhat) / (rhat'*(A(M(rhat))));
        omega = (r'*rhat) / (rhat'*A(rhat));
        x = x + omega * rhat;
        val1 = 0.5*( x'*A(x) ) - x'*b;
        val_record(i+1) = val1;
    end
    %if flag == 0
    %    disp(['richardson iteration failed to converge within ' num2str(maxit) ' iterations.'])
    %end
end
