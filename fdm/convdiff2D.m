function A = convdiff2D(a,b,n,ctype,anisotropic,epsilon)

assert( n > 0 );
assert( ctype == 0 || ctype == 1 || ctype == 2 );
assert( anisotropic == 0 || anisotropic == 1 );
assert( epsilon > 0 );

h = (b-a)/(n+1); 

[y, x] = meshgrid( a+h : h : b-h );
[x, y] = deal( x(:), y(:) );

e = ones(n,1);
I = speye(n,n);
T = spdiags(e*[-1 2 -1],-1:1,n,n) / (h^2); % second derivative
C = spdiags(e*[-1 1],[-1 1],n,n) / (2*h);  % first derivative

% Convection Field
if ( ctype == 1 ) % convection of type 1
    vx = (x-x.^2).*(2*y-1);
    vy = (y-y.^2).*(2*x-1);
end
if ( ctype == 2 ) % convection of type 2
    vx = 4*sin(y).*exp(-x.^2-y.^2).*( cos(x)-2*x.*sin(x) );
    vy = 4*sin(x).*exp(-x.^2-y.^2).*( cos(y)-2*y.*sin(y) );
end

% Diffusion Coefficients (Heterogeneous Diffusion)
dx = zeros(n^2,1);
dy = zeros(n^2,1);

i=1:n^3;

i1 = i( ( x <  1/3 )           & ( y <  1/2 ) );
i2 = i( ( x >= 1/3 & x < 2/3 ) & ( y <  1/2 ) );
i3 = i( ( x >= 2/3 )           & ( y <  1/2 ) );
i4 = i( ( x <  1/3 )           & ( y >= 1/2 ) );
i5 = i( ( x >= 1/3 & x < 2/3 ) & ( y >= 1/2 ) );
i6 = i( ( x >= 2/3 )           & ( y >= 1/2 ) );

dx([i1 i3 i5]) = 1;
dx([i2 i4 i6]) = (anisotropic==1)*1 + (anisotropic==0)*1e3;
dy([i1 i3 i5]) = 1;
dy([i2 i4 i6]) = 1e3;

% matrix assembly
A = dx.*kron(I,T) + dy.*kron(T,I);
A = epsilon*A;

if ( ctype ~= 0 ) % if convection field is of type 1 or 2
    A = A + vx.*kron(I,C) ...
          + vy.*kron(C,I);
end

end