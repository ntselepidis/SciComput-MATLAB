function y = adams_bashforth(f,y1,t,step)

assert((step==3) || (step==4) || (step==5));

n=length(t)-1;
y=zeros(n+1,1);
y(1)=y1;

% runge-kutta (4th order) approximation
for i=1:step-1
    h=t(i+1)-t(i);
    k1=h*f(t(i),y(i));
    k2=h*f(t(i)+0.5*h,y(i)+0.5*k1);
    k3=h*f(t(i)+0.5*h,y(i)+0.5*k2);
    k4=h*f(t(i)+h,y(i)+k3);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
end

if (step==3) 
    % 3-step adams-bashforth
    for i=step:n
        h=t(i+1)-t(i);
        y(i+1)=y(i)+h*(23*f(t(i),y(i))...
                      -16*f(t(i-1),y(i-1))...
                      +5*f(t(i-2),y(i-2)))/12;
    end
elseif (step==4)
    % 4-step adams-bashforth
    for i=step:n
        h=t(i+1)-t(i);
        y(i+1)=y(i)+h*(55*f(t(i),y(i))...
                      -59*f(t(i-1),y(i-1))...
                      +37*f(t(i-2),y(i-2))...
                      -9*f(t(i-3),y(i-3)))/24;
    end
else
    % 5-step adams-bashforth
    for i=step:n
        h=t(i+1)-t(i);
        y(i+1)=y(i)+h*(1901*f(t(i),y(i))...
                      -2774*f(t(i-1),y(i-1))...
                      +2616*f(t(i-2),y(i-2))...
                      -1274*f(t(i-3),y(i-3))...
                      +251*f(t(i-4),y(i-4)))/720;
    end
end

end
