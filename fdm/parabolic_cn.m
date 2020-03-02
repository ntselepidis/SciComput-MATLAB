clear; clc; close all;

alpha   = 1; % diffusion coefficient ( mass / thermal diffusivity )
epsilon = 0; % convection coefficient

n  = 9; % number of unknowns on x
x0 = 0;
x1 = 1;
h = (x1-x0)/(n+1); % mesh size
x = x0+(0:n+1)'*h; % domain inner + boundary points: coordinates

dt = 0.001; % time step
t0 = 0;
t1 = 0.5;
nt = (t1-t0)/dt + 1;

r = (alpha * dt) / h^2;

Afw = spdiags([r/2 1-r r/2].*ones(n,1),(-1:1),n,n);
Abw = spdiags([-r/2 1+r -r/2].*ones(n,1),(-1:1),n,n);

u = zeros(n+2,nt); % matrix holding u at all timesteps (including bounds)
u(:,1)   = sin(pi*x);   % initial condition on all xs
u(1,:)   = u(1,1);      % boundary conditions on x0
u(end,:) = u(1,end);    % boundary conditions on x1
umin = min(u(:,1));
umax = max(u(:,1));
w1 = zeros(n,1); % work vector for boundary conditions (attention: size is n)
w2 = zeros(n,1); % work vector for boundary conditions (attention: size is n)

for t=2:nt
    % compute solution at timestep t
    w1(1) = (r/2)*u(1,  t-1); 
    w1(n) = (r/2)*u(end,t-1);
    w2(1) = (r/2)*u(1,  t);
    w2(n) = (r/2)*u(end,t);
    u(2:end-1,t) = Abw \ (Afw * u(2:end-1,t-1) + w1 + w2); % hybrid scheme
    
    % plot solution at timestep t
    axis([x0 x1 umin umax]); 
    xlabel('x'); ylabel('u');
    getframe; plot(x,u(:,t));
end
