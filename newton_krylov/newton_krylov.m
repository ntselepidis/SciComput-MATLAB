function x = newton_krylov(f,df,xinit,tol,maxit,KSM)
x = xinit;
flag = 0;
PC = KSM.pc; % preconditioner struct
for i=1:maxit

    [A, b] = deal( df(x), f(x) );

    [PC.L, PC.U, PC.P, PC.Q] = factorize( A, PC.All );
    
    dx = ksm(A,b,KSM.tol,KSM.maxit,@(y) KSM.pcfun(y,PC),KSM.type);

    x = x - dx;
    
    if norm(dx) < tol
        flag = 1;
        break;
    end
    
end
if flag 
    disp(['Newton-Krylov method converged in ' num2str(i) ' iterations.']);
else
    disp('Newton-Krylov method failed to converge.')
end
end