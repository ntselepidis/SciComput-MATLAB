% Lu = 0 (2D) on a square domain [0,1]x[0,1]
% neumann on south boundary
% dirichlet on the other boundaries
clear; clc; close all;
u = @(x,y) exp(x).*cos(y);
phi = @(r,c) sqrt(r.^2+c^2);
% r = @(x,xj,y,yj) sqrt((x-xj).^2+(y-yj).^2);
% dux = @(x,xj,y,yj,c) (x-xj) ./ sqrt(r(x,xj,y,yj).^2+c^2);
% du2x = @(x,xj,y,yj,c) ((y-yj).^2+c^2) ./ ((r(x,xj,y,yj,c)+c^2)^(3/2));
c = 1;
a = 0;
b = 1;
nx = 25;
h = 1 / (nx-1);
xo = [a:h:b a:h:b repmat(a, 1, nx-2) repmat(b, 1, nx-2)]';
yo = [repmat(a, 1, nx) repmat(b, 1, nx) a+h:h:b-h a+h:h:b-h]';
[xi, yi] = meshgrid((a+h):h:(b-h)); % inner points
xi = xi(:);
yi = yi(:);
figure, plot(xo, yo, 'r*', xi, yi, '*')
x = [xo; xi];
y = [yo; yi];
n = length(x);
no = length(xo);
r = zeros(n);
dux = zeros(n);
duy = zeros(n);
du2x = zeros(n);
du2y = zeros(n);
for i = 1:n
    j = setdiff(1:n, i);
    r(i,j) = sqrt( (x(i)-x(j)).^2 + (y(i)-y(j)).^2 )';
end
j = 1:n;
for i = 1:n
    dux(i,j) = (x(i)-x(j)) ./ (sqrt(r(i,j).^2+c^2))';
    du2x(i,j) = ((y(i)-y(j)).^2+c^2) ./ ((r(i,j).^2+c^2).^(3/2))';  
    duy(i,j) = (y(i)-y(j)) ./ (sqrt(r(i,j).^2+c^2))';
    du2y(i,j) = ((x(i)-x(j)).^2+c^2) ./ ((r(i,j).^2+c^2).^(3/2))';
end
A = du2x+du2y;
% A(1:no,1:no) = eye(no);
% A(1:no,(no+1):n) = zeros(no,length((no+1):n));
A(1:no,1:n) = phi(r(1:no,1:n),c);
b = zeros(n, 1);
b(1:no) = u(x(1:no),y(1:no));%exp(x(1:no)).*cos(y(1:no));
b(1:nx) = zeros(nx, 1);
A(1:nx,:) = duy(1:nx,:);
figure, imagesc(A)
xx = A \ b;
exact = u(x,y);
xx = phi(r(1:n,1:n), c)*xx;
error = norm(xx-exact) / norm(exact)
figure, plot3(x, y, exact,'*')
figure, plot3(x, y, xx, '*')
