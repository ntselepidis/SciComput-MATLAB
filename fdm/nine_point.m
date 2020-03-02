function A = nine_point(m)

n = (m-1)^2;

e = ones(n,1);

A = spdiags( [e 4*e e 4*e -20*e 4*e e 4*e e], [-m -(m-1) -(m-2) -1 0 1 m-2 m-1 m], n, n );

for i = m-1:m-1:(m-2)*(m-1)
    A(i,i+1) = 0;
    A(i+1,i) = 0;
end

for i = m-1:m-1:(m-3)*(m-1)
    A(i,i+m) = 0;
    A(i+m,i) = 0;
end

for i=1:m-1:(m-1)^2
    A(i,i+m-2) = 0;
    A(i+m-2,i) = 0;
end

end