clear;clc;close all;

f=@(t,y) -20*(y-t^2)+2*t;
y1=-1/3;

h1=0.005;
h2=0.005;
t=[0:h1:0.2 0.2+h2:h2:1]';

y=t.^2-exp(-20*t)/3; % analytical solution

yrk4=runge_kutta_4(f,y1,t);
yab3=adams_bashforth(f,y1,t,3);
yab4=adams_bashforth(f,y1,t,4);
yab5=adams_bashforth(f,y1,t,5);
yam3=adams_moulton(f,y1,t);
ytrp=trapezoid(f,y1,t);

erk4=norm(y-yrk4)/norm(y);
eab3=norm(y-yab3)/norm(y);
eab4=norm(y-yab4)/norm(y);
eab5=norm(y-yab5)/norm(y);
eam3=norm(y-yam3)/norm(y);
etrp=norm(y-ytrp)/norm(y);

figure, 
subplot(2,3,1), 
plot(t,yab3,'ro',t,y,'b'), grid, axis([0 1 -0.5 1]),
title('adams-bashforth-3'), legend(['Error: ' num2str(eab3)]);

% figure,
subplot(2,3,2), 
plot(t,yab4,'ro',t,y,'b'), grid, axis([0 1 -0.5 1]),
title('adams-bashforth-4'), legend(['Error: ' num2str(eab4)]);

% figure,
subplot(2,3,3), 
plot(t,yab5,'ro',t,y,'b'), grid, axis([0 1 -0.5 1]),
title('adams-bashforth-5'), legend(['Error: ' num2str(eab5)]);

% figure,
subplot(2,3,4),
plot(t,yam3,'ro',t,y,'b'), grid, axis([0 1 -0.5 1]),
title('adams-moulton-3'), legend(['Error: ' num2str(eam3)]);

% figure,
subplot(2,3,5),
plot(t,ytrp,'ro',t,y,'b'), grid, axis([0 1 -0.5 1]),
title('trapezoid'), legend(['Error: ' num2str(etrp)]);

% figure,
subplot(2,3,6), 
plot(t,yrk4,'ro',t,y,'b'), grid, axis([0 1 -0.5 1]),
title('runge-kutta-4'), legend(['Error: ' num2str(erk4)]);

