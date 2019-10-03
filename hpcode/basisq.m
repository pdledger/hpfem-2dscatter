function ph=basisq(order,x,y);

% p=0 Edge Functions
ph(1,1) = (2-2.*y)/2.;
ph(1,2) =     0;
  
ph(2,1) =  (-2.*y)/2.;
ph(2,2) =    0;
 
ph(3,1) =   0;
ph(3,2) = (2.*x)/2.;
  
ph(4,1) =   0;
ph(4,2) = (-2+2.*x)/2.;

% p> 0 Edge Functions
if order > 0
ep=[2.*x-1;
   1-2.*x;
   2.*y-1;
   1-2.*y];

dep=[2 0;
     -2 0;
      0 2;
      0 -2];

la=[1-y;
    y;
    x;
    1-x];
dla=[0 -1;
     0 1;
     1 0;
    -1 0];

for p=1:order
for e=1:4
ph(4*p+e,1)=dlegi(ep(e),p-1+2)*dep(e,1)*la(e)+legi(ep(e),p-1+2)*dla(e,1);
ph(4*p+e,2)=dlegi(ep(e),p-1+2)*dep(e,2)*la(e)+legi(ep(e),p-1+2)*dla(e,2);
end
end    

%for p=1:order

%ph(1+4*p,1)=(1-y)*dlegi(2*x-1,p-1+2)*2+legi(2*x-1,p-1+2)*0;

%ph(1+4*p,2)=(1-y)*dlegi(2*x-1,p-1+2)*0+legi(2*x-1,p-1+2)*(-1);

%ph(2+4*p,1)=(y)*dlegi(1-2*x,p-1+2)*(-2)+legi(1-2*x,p-1+2)*0;

%ph(2+4*p,2)=(y)*dlegi(1-2*x,p-1+2)*0+legi(1-2*x,p-1+2)*(1);

%ph(3+4*p,1)=(x)*dlegi(2*y-1,p-1+2)*0+legi(2*y-1,p-1+2)*1;

%ph(3+4*p,2)=(x)*dlegi(2*y-1,p-1+2)*2+legi(2*y-1,p-1+2)*(0);

%ph(4+4*p,1)=(1-x)*dlegi(1-2*y,p-1+2)*0+legi(2*y-1,p-1+2)*(-1);

%ph(4+4*p,2)=(1-x)*dlegi(1-2*y,p-1+2)*(-2)+legi(1-2*y,p-1+2)*(0);

%end
end

edgefun=4*(order+1);

% Interior Functions
if order > 0
basno=edgefun;

% Type 1
for i=0:1:order-1
for j=0:1:order-1
basno=basno+1;
ph(basno,1)=legi(2*y-1,j+2)*dlegi(2*x-1,i+2)*2;
ph(basno,2)=legi(2*x-1,i+2)*dlegi(2*y-1,j+2)*2;
end
end

% Type 2
for i=0:1:order-1
for j=0:1:order-1
basno=basno+1;
ph(basno,1)=legi(2*y-1,j+2)*dlegi(2*x-1,i+2)*1;
ph(basno,2)=-1*legi(2*x-1,i+2)*dlegi(2*y-1,j+2)*1;
end
end

% Type 3
for i=0:1:order-1
basno=basno+1;
ph(basno,1)=legi(2*y-1,i+2);
ph(basno,2)=0;
basno=basno+1;
ph(basno,1)=0;
ph(basno,2)=legi(2*x-1,i+2);
end
end

[m,n]=size(ph);
if m~=4*(order+1)+2*order*(order+1)
  disp('error wrong number of basis functions')
end
  
