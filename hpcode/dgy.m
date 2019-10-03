function out=dgy(x,a,b,x1,x2,y1,y2)

if a> 0 && b> 0
%     this is a curve or an ellipse
%     positive y coordinates
phi1=atan2(y1,x1);
phi2=atan2(y2,x2);

if phi1 < 0
phi1=phi1+2*pi;
end

if phi2 < 0
phi2=phi2+2*pi;
end


%     to account for the fact that 0 and 2pi are the same.      
if abs(phi2-phi1) > pi
   if phi1 > pi
   phi1=phi1-(2*pi);
   else
   phi2=phi2-(2*pi);
   end
end

if abs(phi2-phi1) > pi
disp('error in blending function gy');
end

phi=(phi1*0.5*(1-x))+(phi2*0.5*(1+x));
curve=b*(sin(phi));

u=4.*(curve-(0.5*(1-x)*y1)-(0.5*(1+x)*y2));
du=4.*((b*(cos(phi))*0.5*(phi2-phi1))+(0.5*y1)-(0.5*y2));
v=1-(x^2);
dv=-2*x;

if abs((v^2)) < 0.000001
out=0.;
else
out=((v*du)-(u*dv))/(v^2);
end
 
else
out=0.;
end
  
