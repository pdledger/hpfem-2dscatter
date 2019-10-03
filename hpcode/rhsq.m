function rhsel=rhsq(xyq,order,x,w,rbasisx,rbasisy,lbc,omega,probdata,localblend)

nip=length(x);

esize=4*(order+1)+2*order*(order+1);
rhsel=zeros(esize,1);
bctype=probdata.bctype;

norml=[0 -1; 0 1; 1 0; -1 0]';

for j=1:4
   if lbc(j)~=0 & bctype(lbc(j))==1

% This is a boundary edge and a boundary integral is to be applied   

   for p=1:nip
   
   if j==1
     xi=x(p);
     et=0;
   elseif j==2
     xi=x(p);
     et=1;
   elseif j==3
     xi=1;
     et=x(p);
   else
     xi=0;
     et=x(p);
   end
   phx=rbasisx(:,p+(nip*(j-1)));
   phy=rbasisy(:,p+(nip*(j-1)));
   ph=[phx phy];
   
   J=mapquad(xyq,xi,et,localblend);
   detJ=det(J);
   Jinv=(inv(J))';
   phxy=Jinv*ph';
   phxy=phxy';
   
   nm=inv(J)'*norml(:,j);
   nm=nm./norm(nm);
   
   if j==1 || j==2
   v=sqrt(J(1,1)^2+J(2,1)^2);
   else
   v=sqrt(J(1,2)^2+J(2,2)^2);
   end
   [xp,yp]=getcoord(xyq,xi,et,localblend);

%  wave propagation problems   
%   curl=((kx^2)+(ky^2))*i*(cos((kx*xp)+(ky*yp))-i*sin((kx*xp)+(ky*yp)));
%   curl=-((kx^2)+(ky^2))*exp((kx*xp)+(ky*yp));
%   curl=-exp(1/sqrt(2)*(xp+yp));

    %  use the problem file
    fun=probdata.neufun;
    arg=probdata.neufunarg;
    index=lbc(j);
    tangcurl=fun(xp,yp,index,nm,arg);

   for k=1:esize
     rhsel(k)=rhsel(k)-(phxy(k,1)*tangcurl(1)+phxy(k,2)*tangcurl(2))*w(p)*v;
   end
   end
 end
end
