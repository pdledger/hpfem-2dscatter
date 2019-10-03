function ph=curlbasist(order,x,y);
if order==0
esize=3;
else
esize=(order+1)*(order+2);
end
ph=zeros(esize,1);

% p=0 Edge Functions
ph(1)=2/6*sqrt(3);
ph(2)=2/6*sqrt(3);
ph(3)=2/6*sqrt(3);

% p> 0 Edge Functions
if order > 0
s=[(1/2*y*sqrt(3))-1/2-(1/2*x);
   1/2-(1/2*x)-1/2*y*sqrt(3);
   x];
t=[1/2+1/2*x+1/6*y*sqrt(3);
   1/2-1/2*x+1/6*y*sqrt(3);
   1-1/3*y*sqrt(3)];

ds=[-1/2 1/2*sqrt(3);
    -1/2 -1/2*sqrt(3);
    1 0];

dt=[1/2 1/6*sqrt(3);
    -1/2 1/6*sqrt(3);   
    0 -1/3*sqrt(3)]; 
for i=1:3
  dx(i)=(t(i)*ds(i,1)-s(i)*dt(i,1))/t(i)^2;
  dy(i)=(t(i)*ds(i,2)-s(i)*dt(i,2))/t(i)^2;
end

for p=0:order-1
for e=1:3

%a=-(t(e)^(p-1))*dleg(s(e)/t(e),p-1)*dt(e,1)+(t(e)^(p-1))*dleg(s(e)/t(e),p)*ds(e,1);
%b=-(t(e)^(p-1))*dleg(s(e)/t(e),p-1)*dt(e,2)+(t(e)^(p-1))*dleg(s(e)/t(e),p)*ds(e,2);
%lp=(t(e)^p)*leg(s(e)/t(e),p);

%c=-(t(e)^(p-1+1))*dleg(s(e)/t(e),p-1+1)*dt(e,1)+(t(e)^(p-1+1))*dleg(s(e)/t(e),p+1)*ds(e,1);
%d=-(t(e)^(p-1+1))*dleg(s(e)/t(e),p-1+1)*dt(e,2)+(t(e)^(p-1+1))*dleg(s(e)/t(e),p+1)*ds(e,2);
%lpp=t(e)^(p+1)*leg(s(e)/t(e),p+1);

%ph(3*(p+1)+e)= dt(e,2)*(-a*t(e)-lp*dt(e,1))+...
%               ds(e,2)*c-...
%	       dt(e,1)*(-b*t(e)-lp*dt(e,2))-...
%	       ds(e,1)*d;
ph(3*(p+1)+e)=0;
end
end
%for e=1:3
%for p=1:order
%dse=ds(e,:);
%ph(e+3*p)=ph(e+3*p)/norm(dse);
%end
%end

end
edgefun=3*(order+1);

% Interior Functions
if order > 0
basno=edgefun;
l(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*x-y);
l(2)=y/sqrt(3);
l(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*x-y);
dl=[1/2 -1/(2*sqrt(3));
    0 1/sqrt(3);
    -1/2 -1/(2*sqrt(3))];    
du1x=((l(1)+l(2))*(dl(2,1)-dl(1,1))-(l(2)-l(1))*(dl(1,1)+dl(2,1)))/(l(1)+l(2))^2;
du2x=dl(1,1)+dl(2,1);
du1y=((l(1)+l(2))*(dl(2,2)-dl(1,2))-(l(2)-l(1))*(dl(1,2)+dl(2,2)))/(l(1)+l(2))^2;
du2y=dl(1,2)+dl(2,2);

s=l(2)-l(1);
t=l(1)+l(2);
ds=[dl(2,1)-dl(1,1), dl(2,2)-dl(1,2)];
dt=[dl(2,1)+dl(1,1), dl(2,2)+dl(1,2)];

% Type 1
for i=0:1:order-2
for j=0:1:order-2
if i+j <= order-2
basno=basno+1;

%dL_i^s /dx .../dy
a=-t^(i-1)*dleg(s/t,i-1)*dt(1)+t^(i-1)*dleg(s/t,i)*ds(1);
b=-t^(i-1)*dleg(s/t,i-1)*dt(2)+t^(i-1)*dleg(s/t,i)*ds(2);
lp=t^i*leg(s/t,i);

%dL_i+1^s /dx .../dy
c=-t^(i)*dleg(s/t,i)*dt(1)+t^(i)*dleg(s/t,i+1)*ds(1);
d=-t^(i)*dleg(s/t,i)*dt(2)+t^(i)*dleg(s/t,i+1)*ds(2);
lpp=t^(i+1)*leg(s/t,i+1);
ui=legi(s/t,i+2)*t^(i+2);
%du_i/dx ...dy
duix=-t^(i)*leg(s/t,i)*t*dt(1)+t^(i+1)*leg(s/t,i+1)*ds(1);
duiy=-t^(i)*leg(s/t,i)*t*dt(2)+t^(i+1)*leg(s/t,i+1)*ds(2);

%d/dx(du_i/dy)....d/dy(du_i/dx)
duiyx=dt(2)*(-a*t-lp*dt(1))+ds(2)*(c);
duixy=dt(1)*(-b*t-lp*dt(2))+ds(1)*(d);


vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjxy=dl(3,1)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))...
      +dl(3,2)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1))+...
+l(3)*ddleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))*(dl(3,1)-dl(2,1)-dl(1,1));

dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));
dvjyx=dl(3,2)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1))...
      +dl(3,1)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))+...
+l(3)*ddleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))*(dl(3,1)-dl(2,1)-dl(1,1));


ph(basno)=duiyx*vj+dvjx*duiy+(duix*dvjy+ui*dvjyx)-...
          (duixy*vj+dvjy*duix+(duiy*dvjx+ui*dvjxy));

end
end
end

% Type 2
for i=0:1:order-2
for j=0:1:order-2
if i+j <= order-2
basno=basno+1;

%dL_i^s /dx .../dy
a=-t^(i-1)*dleg(s/t,i-1)*dt(1)+t^(i-1)*dleg(s/t,i)*ds(1);
b=-t^(i-1)*dleg(s/t,i-1)*dt(2)+t^(i-1)*dleg(s/t,i)*ds(2);
lp=t^i*leg(s/t,i);

%dL_i+1^s /dx .../dy
c=-t^(i)*dleg(s/t,i)*dt(1)+t^(i)*dleg(s/t,i+1)*ds(1);
d=-t^(i)*dleg(s/t,i)*dt(2)+t^(i)*dleg(s/t,i+1)*ds(2);
lpp=t^(i+1)*leg(s/t,i+1);
ui=legi(s/t,i+2)*t^(i+2);
%du_i/dx ...dy
duix=-t^(i)*leg(s/t,i)*t*dt(1)+t^(i+1)*leg(s/t,i+1)*ds(1);
duiy=-t^(i)*leg(s/t,i)*t*dt(2)+t^(i+1)*leg(s/t,i+1)*ds(2);

%d/dx(du_i/dy)....d/dy(du_i/dx)
duiyx=dt(2)*(-a*t-lp*dt(1))+ds(2)*(c);
duixy=dt(1)*(-b*t-lp*dt(2))+ds(1)*(d);


vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjxy=dl(3,1)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))...
      +dl(3,2)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1))+...
+l(3)*ddleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))*(dl(3,1)-dl(2,1)-dl(1,1));

dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));
dvjyx=dl(3,2)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1))...
      +dl(3,1)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))+...
+l(3)*ddleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))*(dl(3,1)-dl(2,1)-dl(1,1));


ph(basno)=duiyx*vj+dvjx*duiy-(duix*dvjy+ui*dvjyx)-(...
          duixy*vj+dvjy*duix-(duiy*dvjx+ui*dvjxy));


end
end
end

for j=0:1:order-2
basno=basno+1;
vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));

ph(basno)=2/6*sqrt(3)*vj+...
         1/6*sqrt(3)*(1+x)*dvjx+1/6*y*sqrt(3)*dvjy;
end

end

