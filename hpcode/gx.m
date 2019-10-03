function out=gx(x,a,b,x1,x2,y1,y2)
 
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
disp('error in blending function gx');
end

phi=(phi1*0.5*(1-x))+(phi2*0.5*(1+x));
curve=a*(cos(phi));
      
if abs(1-(x^2)) < 0.000001
out=0.;
else      
out=(4./(1-(x^2)))*(curve-(0.5*(1-x)*x1)-(0.5*(1+x)*x2));
end

else
      
out=0.;
end
