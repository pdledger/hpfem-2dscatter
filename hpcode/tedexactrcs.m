function [rcsval]=tedexactrcs(omega,a,phi,ed,mud,theta)
k = omega; % = omega!
ea = 1;
mua=1;;
kd=sqrt(k^2*ed*mud);
r=a;
rcsval=0;
n=0;
minn=0;
drcsval=1;
while abs(drcsval) > 1e-8 & n <20 | minn==0

    Ja = besselj(n,k*a);%Bessel function of the first kind (a=1)
    Jad = besselj(n,kd*a);
    
    Ya = bessely(n,k*a);%Bessel function of the second kind (a=1)
    
    Ha = Ja - i*Ya;%Hankel function of the second kind
    
    Jr = besselj(n,k*r);%Bessel function of the first kind (a=1)
    
    Jrd = besselj(n,kd*r);
    
    Yr = bessely(n,k*r);%Bessel function of the second kind (a=1)
    
    Hr = Jr - i*Yr;%Hankel function of the second kind
    
    dJa = (besselj(n-1,k*a)-besselj(n+1,k*a))/2;
    
    dJad = (besselj(n-1,kd*a)-besselj(n+1,kd*a))/2;
    
    
    dYa = (bessely(n-1,k*a)-bessely(n+1,k*a))/2;
    
    dHa = dJa - i*dYa;%Hankel function of the second kind
    if n~=0
        F = 2;
    else
        F = 1;
    end
    
    a1 = -Ja/Ha;
    a2 = (mud*dJad)/(mua*kd*a*Jad);
    a3 = -dJa/(k*a*Ja);
    a4 = (mud*dJad)/(mua*kd*a*Jad);
    a5 = -dHa/(k*a*Ha);
    
    
    an = a1*(a2+a3)/(a4+a5);
    
    drcsval=(F*an*cos(n*(phi-theta)));
    rcsval=rcsval+drcsval;
    
    n=n+1;
    if n > 30
    minn=1;
    end
end
    
rcsval=4/k*(abs(rcsval)^2);
