function H = hess2(A)
    [n, m] = size(A);
    if (n ~= m)
        error('Matrix must be square');
    end
    H = A;
    %P = eye(size(A));
    for k = 1:n-2
        x = H(k+1:n, k);
        e1 = eye(length(x), 1);

        if ( x(1) >= 0 )
            sn = 1;
        else
            sn = -1;
        end
        if iscomplex(A)
            sn = exp(j*angle(x(1)));
        end

        u = sn*norm(x)*e1+x;
        u = u / norm(u);

        H(k+1:n, k:n) = H(k+1:n, k:n) - (2*u)*(u'*H(k+1:n, k:n));
        H(1:n, k+1:n) = H(1:n, k+1:n) - (H(1:n, k+1:n)*u)*(2*u)';
        %P(1:n, k+1:n) = P(1:n, k+1:n) - (P(1:n, k+1:n)*u)*(2*u)';
    end
end
