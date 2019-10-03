function [x,w]=gaujac(alf,bet,n)

% MATLAB function for computing points and weights for the integral 
%int_-1^+1 (1-x)^alf(1+x)^bet f(x) dx

maxit=10;

for i=1:n
   if i==1
   an=alf/n;
   bn=bet/n;
   r1=(1+alf)*(2.78/(4+n*n)+0.768*an/n);
   r2=1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
   z=1-r1/r2;
   elseif i==2
   r1=(4.1+alf)/((1+alf)*(1+0.156*alf));
   r2=1+0.06*(n-8)*(1+0.12*alf)/n;
   r3=1+0.12*bet*(1+0.25*abs(alf))/n;
   z=z-(1-z)*r1*r2*r3;
   elseif i==3
   r1=(1.67+0.28*alf)/(1.+0.37*alf);
   r2=1+0.22*(n-8)/n;
   r3=1+8*bet/((6.28+bet)*n*n);
   z=z-(x(1)-z)*r1*r2*r3;
   elseif i==n-1
   r1=(1+0.256*bet)/(0.766+0.119*bet);
   r2=1/(1+0.639*(n-4)/(1+0.71*(n-4)));
   r3=1/(1+0.20*alf/((7.5+alf)*n*n));
   z=z+(z-x(n-3))*r1*r2*r3;
   elseif i==n
   r1=(1+0.37*bet)/(1.67+0.28*bet);
   r2=1/(1+0.22*(n-8)/n);
   r3=1/(1+8*alf/((6.28+alf)*n*n));
   z=z+(z-x(n-2))*r1*r2*r3;
   else
   z=3*x(i-1)-3*x(i-2)+x(i-3);
   end
   alfbet=alf+bet;
   z1=z+0.1;
   its=1;
   while abs(z-z1)>eps && its <=maxit
      temp=2+alfbet;
      p1=(alf-bet+temp*z)/2;
      p2=1;
      for j=2:n
         p3=p2;
	 p2=p1;
	 temp=2*j+alfbet;
	 a=2*j*(j+alfbet)*(temp-2);
	 b=(temp-1)*(alf*alf-bet*bet+temp*(temp-2)*z);
	 c=2*(j-1+alf)*(j-1+bet)*temp;
	 p1=(b*p2-c*p3)/a;
      end
      pp=(n*(alf-bet-temp*z)*p1+2*(n+alf)*(n+bet)*p2)/(temp*(1-z*z));
      z1=z;
      z=z1-p1/pp;
      its=its+1;
   end
   x(i)=z;  
   w(i)=exp(gammaln(alf+n)+gammaln(bet+n)-gammaln(n+1)-gammaln(n+alfbet+1))*temp*2^alfbet/(pp*p2);
end
