function ph=basist(order,x,y);
if order==0
esize=3;
else
esize=(order+1)*(order+2);
end

ph=zeros(esize,2);

% p=0 Edge Functions
ph(1,1)=-1/6*y*sqrt(3);
ph(1,2)=1/6*sqrt(3)*(1+x);

ph(2,1)=-1/6*sqrt(3)*y;
ph(2,2)=1/6*sqrt(3)*(x-1);

ph(3,1)=(1/2)-(1/6*y*sqrt(3));
ph(3,2)=1/6*sqrt(3)*x;

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
%for i=1:3
%  dx(i)=(t(i)*ds(i,1)-s(i)*dt(i,1))/t(i)^2;
%  dy(i)=(t(i)*ds(i,2)-s(i)*dt(i,2))/t(i)^2;
%end

for p=0:order-1
for e=1:3

if t(e)==0
 a=0; %the limit of ls(s,t,n) =0 for n>=2 but a=t(e)*ls(s,t,n) for n>=0 and so a=0 for all n when t(e)=0
 if p>=1
   b=0;  % limit of ls(s,t,n) =0 for n >=2
 else
   b=s(e)/2 % limit only true for n>=2 but for n=1 ls(s,t,1)=s
 end
else
 a=t(e)^(p+1)*leg(s(e)/t(e),p);
 b=t(e)^(p+1)*leg(s(e)/t(e),p+1);
end
ph(3*(p+1)+e,1)=-a*dt(e,1)+b*ds(e,1);
ph(3*(p+1)+e,2)=-a*dt(e,2)+b*ds(e,2);
%p*t(e)^(p-1)*dt(e,1)*legi(s(e)/t(e),p)+t(e)^p*dlegi(s(e)/t(e),p)*dx(e);
%ph(3*(p+1)+e,2)=p*t(e)^(p-1)*dt(e,2)*legi(s(e)/t(e),p)+t(e)^p*dlegi(s(e)/t(e),p)*dy(e);
end
end    


%for e=1:3
%for p=1:order
%dse=ds(e,:);
%ph(e+3*p,1)=ph(e+3*p,1)/norm(dse);
%ph(e+3*p,2)=ph(e+3*p,2)/norm(dse);
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

s=l(2)-l(1);
t=l(1)+l(2);
ds=[dl(2,1)-dl(1,1), dl(2,2)-dl(1,2)];
dt=[dl(2,1)+dl(1,1), dl(2,2)+dl(1,2)];

% Type 1
for i=0:1:order-2
for j=0:1:order-2
if i+j <= order-2
basno=basno+1;

if t==0
  if i>= 1
    duix=0;
    duiy=0;
    ui=0;
  else
    duix=s*ds(1);
    duiy=s*ds(2);
    ui=s^2/2;
  end
else
  ui=legi(s/t,i+2)*t^(i+2);
  duix=-leg(s/t,i)*t^i*t*dt(1)+leg(s/t,i+1)*t^(i+1)*ds(1);
  duiy=-leg(s/t,i)*t^i*t*dt(2)+leg(s/t,i+1)*t^(i+1)*ds(2);
end

vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));

ph(basno,1)=vj*duix+ui*dvjx;
ph(basno,2)=vj*duiy+ui*dvjy;
end
end
end

% Type 2
for i=0:1:order-2
for j=0:1:order-2
if i+j <= order-2
basno=basno+1;
if t==0
  if i>= 1
    duix=0;
    duiy=0;
    ui=0;
  else
    duix=s*ds(1);
    duiy=s*ds(2);
    ui=s^2/2;
  end
else
  ui=legi(s/t,i+2)*t^(i+2);
  duix=-leg(s/t,i)*t^i*t*dt(1)+leg(s/t,i+1)*t^(i+1)*ds(1);
  duiy=-leg(s/t,i)*t^i*t*dt(2)+leg(s/t,i+1)*t^(i+1)*ds(2);
end

vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));

ph(basno,1)=vj*duix-ui*dvjx;
ph(basno,2)=vj*duiy-ui*dvjy;
end
end
end

% Type 3
for j=0:1:order-2

basno=basno+1;
vj=l(3)*leg(l(3)-l(2)-l(1),j);
ph(basno,1)=-1/6*y*sqrt(3)*vj;
ph(basno,2)=1/6*sqrt(3)*(1+x)*vj;
end
end

[m,n]=size(ph);
if m~=esize
  disp('error wrong number of basis functions')
end

