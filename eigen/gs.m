function Q=gs(A)
    n=size(A,1);
    k=size(A,2);
    Q=zeros(n,k);
    Q(:,1)=A(:,1)/sqrt(A(:,1)'*A(:,1));
    for i=2:k
        Q(:,i)=A(:,i);
        for j=1:i-1
            Q(:,i)=Q(:,i)-( Q(:,i)'*Q(:,j) )/( Q(:,j)'*Q(:,j) )*Q(:,j);
        end
        Q(:,i)=Q(:,i)/sqrt(Q(:,i)'*Q(:,i));
    end
end