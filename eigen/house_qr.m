function [Q,R]=house_qr(A)
[m,n]=size(A);
Q=eye(m);
R=A;
for k=1:n-1
    x=R(k:m,k);
    e=eye(length(x),1);
    
    if (x(1)>=0)
        sn=1;
    else
        sn=-1;
    end
    
    u=sn*norm(x)*e+x;
    u=u/norm(u);
    
    R(k:m,:) = R(k:m,:) - (2*u)*(u'*R(k:m,:));
    Q(:,k:m) = Q(:,k:m) - (Q(:,k:m)*u)*(2*u)';
end
end