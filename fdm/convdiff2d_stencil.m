clear; clc; close all;

% Spatial domain data
[xmin, xmax, nx] = deal(0.0, 1.0, 100);
[ymin, ymax, ny] = deal(0.0, 1.0, 100);

hx = (xmax - xmin) / (nx - 1);
hy = (ymax - ymin) / (ny - 1);

% Stream function and velocities vx, vy
B = 100;

psi = @(x, y) B * sin((pi * x) / xmax) .* sin((pi * y) / ymax);

[y, x] = meshgrid(ymin:hy:ymax, xmin:hx:xmax);

S = psi(x, y);

% Sx
Sx = zeros(nx, ny);
% Sx(1, :) = (S(2, :) - S(1, :)) / hx;
% Sx(nx, :) = (S(nx, :) - S(nx-1, :)) / hx;
Sx(2:nx-1, :) = (S(3:nx, :) - S(1:nx-2, :)) / (2 * hx);

% Sy
Sy = zeros(nx, ny);
% Sy(:, 1) = (S(:, 2) - S(:, 1)) / hy;
% Sy(:, nx) = (S(:, ny) - S(:, ny-1)) / hy;
Sy(:, 2:ny-1) = (S(:, 3:ny) - S(:, 1:ny-2)) / (2 * hy);

% vx, vy
[vx, vy] = deal( Sy, -Sx );

figure, contourf(x, y, S), hold on, quiver(x, y, Sx, Sy), hold off;
figure, contourf(x, y, S), hold on, quiver(x, y, vx, vy), hold off;

vxmax = max(vx, [], 'all');
vymax = max(vy, [], 'all');

% Other params, CFL, and timestep configuration

k = 1.0;
a_dif = 0.2;
a_adv = 0.4;
total_time = 10.0;

dt_dif = ( a_dif * min(hx, hy)^2 ) / k
dt_adv = a_adv * min(hx / vxmax, hy / vymax)
dt = min( dt_dif, dt_adv )
nsteps = ceil(total_time / dt)

% Initial state
T = zeros(nx, ny);
dT2 = zeros(nx, ny);
dTx = zeros(nx, ny);
dTy = zeros(nx, ny);

figure,
for step = 1 : nsteps
    % Dirichlet BCs on bottom and top boundaries
    T(:, 1) = 1.0;
    T(:, ny) = 0.0;
    % Neumann BCs on left and right boundaries
    T(1, :) = T(2, :);
    T(nx, :) = T(nx-1, :);

    % Diffusion term
    dT2(2:nx-1, 2:ny-1) = k * ( T(3:nx,   2:ny-1) ...
                              + T(1:nx-2, 2:ny-1) ...
                              + T(2:nx-1, 3:ny  ) ...
                              + T(2:nx-1, 1:ny-2) ...
                              - 4.0 * T(2:nx-1, 2:ny-1) ) / hx^2;
    % Advection term
    posx  = vx(2:nx-1, 2:ny-1) > 0;
    dTxBW = vx(2:nx-1, 2:ny-1) .* (T(2:nx-1, 2:ny-1) - T(1:nx-2, 2:ny-1)) / hx;
    dTxFW = vx(2:nx-1, 2:ny-1) .* (T(3:nx,   2:ny-1) - T(2:nx-1, 2:ny-1)) / hx;
    dTx(2:nx-1, 2:ny-1) = posx .* dTxBW + (~posx) .* dTxFW;

    posy  = vy(2:nx-1, 2:ny-1) > 0;
    dTyBW = vy(2:nx-1, 2:ny-1) .* (T(2:nx-1, 2:ny-1) - T(2:nx-1, 1:ny-2)) / hy;
    dTyFW = vy(2:nx-1, 2:ny-1) .* (T(2:nx-1, 3:ny  ) - T(2:nx-1, 2:ny-1)) / hy;
    dTy(2:nx-1, 2:ny-1) = posy .* dTyBW + (~posy) .* dTyFW;

    % Integration
    T = T + dt * ( dT2 - dTx - dTy );

    % Plot solution
    getframe; contourf(x, y, T), title(strcat("step ", num2str(step)))

end
