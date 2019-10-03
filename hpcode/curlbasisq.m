function ph=curlbasisq(order,x,y);

ph=zeros(4*(order+1)+2*order*(order+1),1);
% p=0 Edge Functions
ph(1)=1;
ph(2)=1;
ph(3)=1;
ph(4)=1;

% p> 0 Edge Functions
if order > 0
ep=[2*x-1,  1-2*x,   2*y-1,   1-2*y];

dep=[2 0;
     -2 0;
      0 2;
      0 -2];

la=[1-y,    y,    x,    1-x];
dla=[0 -1;
     0 1;
     1 0;
    -1 0];

for p=1:order
for e=1:4

%a=dep(e,2).*(ddlegi(ep(e),p-1+2).*dep(e,1).*la(e)+dlegi(ep(e),p-1+2).*dla(e,1))+...
%  dla(e,2).*(dlegi(ep(e),p-1+2).*dep(e,1));
%b=dep(e,1).*(ddlegi(ep(e),p-1+2).*dep(e,2).*la(e)+dlegi(ep(e),p-1+2).*dla(e,2))+...
%  dla(e,1).*(dlegi(ep(e),p-1+2).*dep(e,2));  

ph(e+4*p)=0;%a-b;
end

end
end

edgefun=4*(order+1);

% Interior Functions
if order > 0
basno=edgefun;

% Type 1
for i=0:1:order-1
for j=0:1:order-1
basno=basno+1;
ph(basno)=0;
end
end

% Type 2
for i=0:1:order-1
for j=0:1:order-1
basno=basno+1;
ph(basno)=-4*dlegi(2*x-1,i+2).*dlegi(2*y-1,j+2);
end
end

% Type 3
for i=0:1:order-1
basno=basno+1;
ph(basno)=-2*dlegi(2*y-1,i+2);
basno=basno+1;
ph(basno)=2*dlegi(2*x-1,i+2);
end
end

m=length(ph);
if m~=4*(order+1)+2*order*(order+1)
  disp('error wrong number of curl basis functions')
end


