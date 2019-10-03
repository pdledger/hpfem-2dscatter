function [dirval]=dirichlett(edges,xt,wt,glob,coord,intmat,unkl,unkh,dir,...                        
		rbasisxt,rbasisyt,nelemt,nelemq,dirval,order,probdata,coefft)

nip=length(xt);
bctype=probdata.bctype;

esize=3;
if order > 0
esize=(order+1)*(order+2);
end

norml=[0.5*sqrt(3.) 0.5; -0.5*sqrt(3.) 0.5; 0 -1.]';
for i=1:nelemt

   for l=1:3
      for k=1:2
         xyt(l,k)=coord(intmat(i,l),k);
      end
   end
   
   for l=1:3
      for k=1:2
      localblendt(l,k)=coefft(i,l,k);
      end
   end


   for jj=1:3
   if edges(glob(i+nelemq,jj),3)~=0 & bctype(edges(glob(i+nelemq,jj),3))==2
%  then we have a Dirichlet edge
   a=zeros(order+1,order+1);
   r=zeros(order+1,1);

   for p=1:nip
    if jj==1
     eta=sqrt(3.)*0.5*(xt(p)+1);
     xi=(-1+xt(p))*(-0.5);
   elseif jj==2
     eta=sqrt(3.)*(-0.5)*(xt(p)-1);
     xi=-1*(1+xt(p))*0.5;
   else
    eta=0.;
    xi=xt(p);
   end
   phx=rbasisxt(:,p+(nip*(jj-1)));
   phy=rbasisyt(:,p+(nip*(jj-1)));
   ph=[phx phy];

 
   J=maptri(xyt,xi,eta,localblendt);
   detJ=det(J);
   Jinv=(inv(J))';
   phxy=Jinv*ph';
   phxy=phxy';      
   
   nm=inv(J)'*norml(:,jj);
   nm=nm./norm(nm);
   
   if jj==1
     ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*0.5));
     ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*0.5));
     v=sqrt((ddgx^2)+(ddgy^2));
   elseif jj==2
      ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*-0.5));
      ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*-0.5));
      v=sqrt((ddgx^2)+(ddgy^2));
   else
      ddgx=J(1,1);
      ddgy=J(2,1);
      v=sqrt((ddgx^2)+(ddgy^2));
   end 
   [xp,yp]=getcoordt(xyt,xi,eta,localblendt);


    %  use the problem file
    fun=probdata.dirfun;
    arg=probdata.dirfunarg;
    index=edges(glob(i+nelemq,jj),3);
    tange=fun(xp,yp,index,nm,arg);


%  the tangential components of all basis functions except those
%  on edge j vanish   
   for q=1:order+1
%   r(q)=r(q)+(wt(p)*v*(nm(1)*e(2)-nm(2)*e(1))*(nm(1)*phxy(jj+(3*(q-1)),2)-nm(2)*phxy(jj+(3*(q-1)),1)));
%   For Scattering Problems n x E^s = - n x E^i
   r(q)=r(q)+wt(p)*v*tange*(nm(1)*phxy(jj+(3*(q-1)),2)-nm(2)*phxy(jj+(3*(q-1)),1));; 
   end
   for q=1:order+1
   for qq=1:order+1
   a(q,qq)=a(q,qq)+(wt(p)*v*(nm(1)*phxy(jj+(3*(q-1)),2)-nm(2)*phxy(jj+(3*(q-1)),1))*...
   (nm(1)*phxy(jj+(3*(qq-1)),2)-nm(2)*phxy(jj+(3*(qq-1)),1)));
   end
   end
   end
   if dir(i+nelemq,jj) < 0
     disp('wrong orientation defined in mesh!')
   end
   if det(a)==0
   disp('error cannot apply Boundary conditions')
   end
   s=a\r;
   % lowest order block
   dirval(abs(unkl(glob(i+nelemq,jj))))=s(1)*dir(i+nelemq,jj);
   %  higher order block
   for q=1:order
   dirval(abs(unkh(glob(i+nelemq,jj),q)))=s(q+1)*(dir(i,jj)^(q-1));
   end
   end
   end
end   
