function [x,y]=getcoordt(xy,xi,eta,localblend)
l(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*xi-eta);
l(2)=eta/sqrt(3);
l(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*xi-eta);
x=0;
y=0;
for j=1:3
  x=x+(xy(j,1)*l(j));
  y=y+(xy(j,2)*l(j));
end

      x1=(l(1)*l(2)*gx(l(2)-l(1),localblend(1,1),localblend(1,2)...
,xy(1,1),xy(2,1),xy(1,2),xy(2,2)));
      y1=(l(1)*l(2)*gy(l(2)-l(1),localblend(1,1),localblend(1,2)...
,xy(1,1),xy(2,1),xy(1,2),xy(2,2)));
      
      
      x2=(l(2)*l(3)*gx(l(3)-l(2),localblend(2,1),localblend(2,2)...
,xy(2,1),xy(3,1),xy(2,2),xy(3,2)));
      y2=(l(2)*l(3)*gy(l(3)-l(2),localblend(2,1),localblend(2,2)...
,xy(2,1),xy(3,1),xy(2,2),xy(3,2)));
      

      x3=(l(3)*l(1)*gx(l(1)-l(3),localblend(3,1),localblend(3,2)...
,xy(3,1),xy(1,1),xy(3,2),xy(1,2)));
      y3=(l(3)*l(1)*gy(l(1)-l(3),localblend(3,1),localblend(3,2)...
,xy(3,1),xy(1,1),xy(3,2),xy(1,2)));


      x=x+x1+x2+x3;
      y=y+y1+y2+y3;
