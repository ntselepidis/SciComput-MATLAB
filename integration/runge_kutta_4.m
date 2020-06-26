function y = runge_kutta_4(f,y1,t)

n=length(t)-1;
y=zeros(n+1,1);
y(1)=y1;

% runge-kutta (4th order) approximation
for i=1:n
    h=t(i+1)-t(i);
    k1=h*f(t(i),y(i));
    k2=h*f(t(i)+0.5*h,y(i)+0.5*k1);
    k3=h*f(t(i)+0.5*h,y(i)+0.5*k2);
    k4=h*f(t(i)+h,y(i)+k3);
    y(i+1)=y(i)+(k1+2*k2+2*k3+k4)/6;
end

end