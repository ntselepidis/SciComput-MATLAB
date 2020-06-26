function y = trapezoid(f,y1,t)

n=length(t)-1;
y=zeros(n+1,1);
y(1)=y1;

for i=1:n
    h=t(i+1)-t(i);
    y(i+1)=y(i)+h*f(t(i),y(i)); % euler predictor
    y(i+1)=y(i)+h*(f(t(i),y(i))+f(t(i+1),y(i+1)))/2; % trapezoid corrector
end

end