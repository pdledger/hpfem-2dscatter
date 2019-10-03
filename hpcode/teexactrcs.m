function [rcsval]=teexactrcs(omega,a,phi,theta)
x=omega*a;
rcsval=0;
n=0;
drcsval=1;
minn=0;
minnval=30;
while abs(drcsval) > 1e-8 | minn==0
dj=(-besselj(n-1,x)+besselj(n+1,x))/2;
dy=(-bessely(n-1,x)+bessely(n+1,x))/2;
h=dj-sqrt(-1)*dy;
if n==0
drcsval=(dj/h)*cos(n*phi);
else
drcsval=2*(dj/h)*cos(n*(phi-theta));
end
rcsval=rcsval+drcsval;
n=n+1;
if n > minnval
minn=1;
end
end
rcsval=4/omega*(abs(rcsval)^2);

