function [rcsval]=tmexactrcs(omega,a,phi,theta)
x=omega*a;
rcsval=0;
n=0;
minn=0;
drcsval=1;
minnval=30;
while abs(drcsval) > 1e-8 | minn==0
j=besselj(n,x);
y=bessely(n,x);
h=j-sqrt(-1)*y;
if n==0
drcsval=(j/h)*cos(n*phi);
else
drcsval=2*(j/h)*cos(n*(phi-theta));
end
rcsval=rcsval+drcsval;
n=n+1;
if n > minnval
minn=1;
end
end
rcsval=4/omega*(abs(rcsval)^2);
