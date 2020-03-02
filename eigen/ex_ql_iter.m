clear;clc;format long;

A = [8 0.25 0.5 2 -1; 
     0.25 -4 0 1 2; 
     0.5 0 5 0.75 -1; 
     2 1 0.75 5 -0.5; 
     -1 2 -1 -0.5 6];

flag = 1;

[vec_approx,val_approx]=ql_iter(A,1e-10,200);

if (flag == 1)
    disp(' ');
    [vec,val]=eig(A);
    [eigval1,p1]=sort(diag(val));
    [eigval2,p2]=sort(val_approx);
    disp('Exact and approximate eigenvalues:'); 
    disp([eigval1 eigval2]);
    disp('Exact eigenvectors:');
    disp(vec(:,p1));
    disp('Approximate eigenvectors:');
    disp(vec_approx(:,p2));
end

