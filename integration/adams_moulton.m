function y = adams_moulton(f,y1,t)

n=length(t)-1;
y=zeros(n+1,1);
y(1)=y1;
step=4;

% runge-kutta (4th order) approximation
for i=1:step-1
    h=t(i+1)-t(i);
    k1=h*f(t(i),y(i));
    k2=h*f(t(i)+0.5*h,y(i)+0.5*k1);
    k3=h*f(t(i)+0.5*h,y(i)+0.5*k2);
    k4=h*f(t(i)+h,y(i)+k3);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
end

for i=step:n
    h=t(i+1)-t(i);
    % 4-step adams-bashforth
    y(i+1)=y(i)+h*(55*f(t(i),y(i))...
                  -59*f(t(i-1),y(i-1))...
                  +37*f(t(i-2),y(i-2))...
                  -9*f(t(i-3),y(i-3)))/24;
    % 3-step adams-moulton
    y(i+1)=y(i)+h*(9*f(t(i+1),y(i+1))...
                 +19*f(t(i),y(i))...
                  -5*f(t(i-1),y(i-1))...
                    +f(t(i-2),y(i-2)))/24;
end


end