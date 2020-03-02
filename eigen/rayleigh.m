function [x,q]=rayleigh(A,x,q,tol,maxit)
I = speye(length(A));
for i = 1:maxit
    x = ( A - q*I )\x;
    x = x/norm(x);
    y = A*x;
    q = x'*y; % rayleigh quotient: (x'*A*x)/(x'*x)
    r = y - q*x;
    if ( norm(r) < tol )
        fprintf("Rayleigh iteration converged in %d iterations.\n",i);
        break;
    end
end
end