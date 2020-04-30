% Lu = 0 (2D) on a square domain [0,1]x[0,1]
% neumann on south boundary
% dirichlet on the other boundaries
clear; clc; close all;
rng(0); % used for reproducibility
u = @(x,y) exp(x).*cos(y);
phi = @(r,c) sqrt(r.^2+c^2);
c = 1; % shape parameter
a = 0;
b = 1;
nx = 25;
h = 1 / (nx-1);
xo = [a:h:b a:h:b repmat(a, 1, nx-2) repmat(b, 1, nx-2)]';
yo = [repmat(a, 1, nx) repmat(b, 1, nx) a+h:h:b-h a+h:h:b-h]';
% [xi, yi] = meshgrid((a+h):h:(b-h)); % inner points
% xi = xi(:);
% yi = yi(:);
xiyi = rand((nx-4)*nx, 2);
xiyi(xiyi(:,1)==a,:) = [];
xiyi(xiyi(:,1)==b,:) = [];
xiyi(xiyi(:,2)==a,:) = [];
xiyi(xiyi(:,2)==b,:) = [];
xi = xiyi(:,1); yi = xiyi(:,2);
figure, plot(xo, yo, 'r*', xi, yi, '*'), grid
x = [xo; xi];
y = [yo; yi];
n = length(x);
no = length(xo);
r = zeros(n);
dphix = zeros(n);
dphiy = zeros(n);
dphi2x = zeros(n);
dphi2y = zeros(n);
for i = 1:n
    j = setdiff(1:n, i);
    r(i,j) = sqrt( (x(i)-x(j)).^2 + (y(i)-y(j)).^2 )';
end
j = 1:n;
for i = 1:n
    dphix(i,j) = (x(i)-x(j)) ./ (sqrt(r(i,j).^2+c^2))';
    dphiy(i,j) = (y(i)-y(j)) ./ (sqrt(r(i,j).^2+c^2))';
    dphi2x(i,j) = ((y(i)-y(j)).^2+c^2) ./ ((r(i,j).^2+c^2).^(3/2))';  
    dphi2y(i,j) = ((x(i)-x(j)).^2+c^2) ./ ((r(i,j).^2+c^2).^(3/2))';
end
A = dphi2x+dphi2y;
b = zeros(n, 1);
% Dirichlet Boundaries
A(1:no,:) = phi(r(1:no,:),c);
b(1:no) = u(x(1:no),y(1:no));
% Neumann Boundaries (overwrite south boundary)
A(1:nx,:) = dphiy(1:nx,:);
b(1:nx) = zeros(nx, 1);
% Solve
w = A \ b;
xx = phi(r,c)*w;
exact = u(x,y);
error = norm(xx-exact) / norm(exact)
% figure, imagesc(A);
figure, plot3(x, y, exact, '*'), grid, title('Exact solution');
figure, plot3(x, y, xx, '*'), grid, title('Meshless solution');
% figure, plot3(x(1:no),y(1:no),xx(1:no),'r*', ...
%               x(no+1:end),y(no+1:end),xx(no+1:end),'b*')
