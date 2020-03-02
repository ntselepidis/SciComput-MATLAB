function A = five_point(m)

n = (m-1)^2;

e = ones(n,1);

A = spdiags( [e e -4*e e e], [-(m-1) -1 0 1 m-1], n, n );

for i=1:m-2
    A(i*(m-1),i*(m-1)+1) = 0;
    A(i*(m-1)+1,i*(m-1)) = 0;
end

end