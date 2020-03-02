clear; clc; close all;

alpha   = 1; % diffusion coefficient ( mass / thermal diffusivity )
epsilon = 0; % convection coefficient
fw      = 0; % 1 for fw diff (explicit) || 0 for bw diff (implicit)

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
I = speye(n,n);
D = spdiags([1 -2 1].*ones(n,1),(-1:1),n,n);
C = spdiags([-1 1]  .*ones(n,1),[-1 1],n,n);

if (fw)
%     assert( r <= 0.5 ); % check for stability of explicit method
%     A = spdiags([r 1-2*r r].*ones(n,1),(-1:1),n,n);
    A = I + D*(alpha*dt)/h^2 + C*(epsilon*dt)/(2*h);
else
%     A = spdiags([-r 1+2*r -r].*ones(n,1),(-1:1),n,n);
    A = I - D*(alpha*dt)/h^2 - C*(epsilon*dt)/(2*h);
end

u = zeros(n+2,nt); % matrix holding u at all timesteps (including bounds)
u(:,1)   = sin(pi*x);   % initial condition on all xs
u(1,:)   = u(1,1);      % boundary conditions on x0
u(end,:) = u(1,end);    % boundary conditions on x1
umin = min(u(:,1));
umax = max(u(:,1));
w = zeros(n,1); % work vector for boundary conditions (attention: size is n)

for t=2:nt
    % compute solution at timestep t
    if (fw)
        w(1) = r*u(1,  t-1);
        w(n) = r*u(end,t-1);
        u(2:end-1,t) = A * u(2:end-1,t-1) + w;
    else
        w(1) = r*u(1,  t);
        w(n) = r*u(end,t);
        u(2:end-1,t) = A \ ( u(2:end-1,t-1) + w );
    end
    
    % plot solution at timestep t
    axis([x0 x1 umin umax]); 
    xlabel('x'); ylabel('u');
    getframe; plot(x,u(:,t));
end
