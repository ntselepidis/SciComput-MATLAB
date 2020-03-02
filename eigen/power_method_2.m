function [x,lambda] = power_method_2(A,x,lambda,tol,maxit,nrm_t)
if (nargin==5)
    nrm_t=2;
end

y = feval(A,x);
lambda = norm(y,nrm_t);
for i=1:maxit
    r = y-lambda*x;
    if (norm(r,nrm_t)/norm(y,nrm_t) < tol)
        fprintf("Power method converged in %d iterations.\n",i);
        break;
    end
    x = y/lambda;
    y = feval(A,x);
    lambda = norm(y,nrm_t);
end
end