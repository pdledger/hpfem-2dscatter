function [epett,layer]=splittri(order2,nptri)

layer=zeros(nptri+1,nptri);
epett=zeros(nptri,2);

nlayer=order2+2;
dy=sqrt(3.)/real(nlayer);      
np=0;
for i=0:nlayer
% calculate y for layer
   y=sqrt(3.)-(i*dy);
   if i==0
      dx=0.;
   else
      dx=(-(2*y)+(2*sqrt(3.)))/((sqrt(3.))*i);
   end
   for j=1:1:i+1
% calculate x coord
      x=((y-sqrt(3.))/sqrt(3.))+(dx*(j-1));
      np=np+1;
      layer(i+1,j)=np;
      epett(np,1)=x;
      epett(np,2)=y;
   end
end


