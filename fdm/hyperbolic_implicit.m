clear; clc; close all;

L = 1;
n = 100;
h = L/(n+1);

T = 4;
nT = 500;
dt = T/nT;

r = (dt/h)^2;
% assert(r<=1);

A = spdiags(ones(n,1)*[-r 4+2*r -r], -1:1, n, n);
B = spdiags(ones(n,1)*[2*r 4*(2-r) 2*r], -1:1, n, n);

x = (0:h:L)';
u = zeros(n+2,nT);
u(:,1) = sin(2*pi*x)*0;
% v = 2*pi*sin(2*pi*x);
v = (1-cos(2*pi*x))/2;
% m = max(u(:,1))*2;
m = 0.5;

u(2:end-1,2) = (2*A) \ ( B*u(2:end-1,1) + A*(2*dt*v(2:end-1)) );
for t = 3:nT
    u(2:end-1,t) = A \ ( B*u(2:end-1,t-1) - A*u(2:end-1,t-2) );
    plot(x, u(:,t));
    axis([0 L -m m]); 
    grid;
    getframe;
end

% figure, mesh(u);
