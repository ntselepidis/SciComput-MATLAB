clear;clc;format long;

A = [3 1 1; -1 3 1];

flag = 1;

[U,SU]=qr_iter(A*A',1e-16,150);
[V,SV]=qr_iter(A'*A,1e-16,150);

S=U'*A*V;
sn=sign( diag(S) )'; % resolving sign ambiguity: singular values > 0
S=abs(S);

if (flag == 1)
    disp(' ');
    [UU,SS,VV]=svd(A);
    disp(UU*SS*VV');
%     disp(U*S*V');
    disp((sn.*U)*S*V');
    norm(A-(sn.*U)*S*V')
end

