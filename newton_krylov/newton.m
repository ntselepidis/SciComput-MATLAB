function x = newton(f,df,xinit,tol,maxit)
x = xinit;
flag = 0;
for i=1:maxit
    dx = df(x) \ f(x);
    x = x - dx;
    if norm(dx) < tol
        flag = 1;
        break;
    end
end
if flag 
    disp(['Newton method converged in ' num2str(i) ' iterations.']);
%     disp('x = '); disp(x);
else
    disp('Newton method failed to converge.')
end
end