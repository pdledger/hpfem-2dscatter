function [x,w]=gaussquad(x1,x2,n)      
% Compute the weights and locations for 1--dimensional gauss quadrature
m=(n+1)/2;
xm=0.5*(x2+x1);;
xl=0.5*(x2-x1);
for i=1:m
z=cos(pi*(i-0.250)/(n+0.50));
z1=z+0.1;
while abs(z-z1)>eps
p1=1.0;
p2=0.0;
for j=1:n
p3=p2;
p2=p1;
p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
end
pp=n*(z*p1-p2)/(z*z-1.0);
z1=z;
z=z1-p1/pp;
end
x(i)=xm-xl*z;
x(n+1-i)=xm+xl*z;
w(i)=2.d0*xl/((1.0-z*z)*pp*pp);
w(n+1-i)=w(i);
end
