function [x,lambda] = power_method(A,x,lambda,tol,maxit,nrm_t)
if (nargin==5)
    nrm_t=2;
end
for i=1:maxit
    y = feval(A,x);
    lambda_old = lambda;
    lambda = norm(y,nrm_t);
    x = y/lambda;
    r = lambda - lambda_old;
    if (abs(r)/abs(lambda) < tol)
        fprintf("Power method converged in %d iterations.\n",i);
        break;
    end
end
end