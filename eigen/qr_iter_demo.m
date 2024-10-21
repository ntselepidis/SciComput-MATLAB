clear;clc;format long;

A = [8 0.25 0.5 2 -1; 
     0.25 -4 0 1 2; 
     0.5 0 5 0.75 -1; 
     2 1 0.75 5 -0.5; 
     -1 2 -1 -0.5 6];

flag = 1;
maxit = 500;
tol = 1e-10;
hessenberg_reduction = 1;
single_shift = 1;

if (hessenberg_reduction == 0)
    disp('Skipping Hessenberg reduction ...')
    [vec_approx,val_approx]=qr_iter(A,tol,maxit);
else
    disp('Performing Hessenberg reduction ...')
    [P, H] = hess(A);
    if (single_shift == 0)
        [vec_approx,val_approx]=qr_iter2(H,tol,maxit);
    else
        [vec_approx,val_approx]=qr_iter2_with_rayleigh_quotient_shift(H,tol,maxit);
    end
end

if (flag == 1)
    disp(' ');
    [vec,val]=eig(A);
    [eigval1,p1]=sort(diag(val));
    [eigval2,p2]=sort(val_approx);
    disp('Exact and approximate eigenvalues:'); 
    disp([eigval1 eigval2]);
%    disp('Exact eigenvectors:');
%    disp(vec(:,p1));
%    disp('Approximate eigenvectors:');
%    disp(vec_approx(:,p2));
end

%disp('------------------');
%disp('Reconstructing ...');
%disp('------------------');
%A
%if (flag == 1)
%    A_tilde1 = vec*val*vec'
%end
%if (hessenberg_reduction == 0)
%    A_tilde2 = vec_approx*diag(val_approx)*vec_approx'
%else
%    A_tilde2 = P*(vec_approx*diag(val_approx)*vec_approx')*P'
%    %A_tilde2 = (P*vec_approx)*diag(val_approx)*(P*vec_approx)'
%end
%
%disp('---------------------');
%disp('Matrix difference ...');
%disp('---------------------');
%if (flag == 1)
%    A - A_tilde1
%end
%A - A_tilde2
