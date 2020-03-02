function x = ksm(A,b,tol,maxit,M,krylov_sel)

if ( krylov_sel == 0 )
    x = bicgstab(A,b,tol,maxit,M);
elseif ( krylov_sel == 1 )
    x = gmres(A,b,20,tol,maxit,M);
else
    x = gmres(A,b,[],tol,maxit,M);
end

end