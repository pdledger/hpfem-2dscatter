function rhsel=rhst(xyt,order,xt,wt,rbasisxt,rbasisyt,lbc,probdata,localblendt,epr,mur,...
curlbast,basisxt,basisyt,intw,intxi,inteta,omega,lfnum)

nip=length(xt);
if order==0
esize=3;
else
esize=(order+1)*(order+2);
end

rhsel=zeros(esize,1);
bctype=probdata.bctype;

norml=[0.5*sqrt(3.) 0.5; -0.5*sqrt(3.) 0.5; 0 -1.]';
for j=1:3
   if lbc(j)~=0 & bctype(lbc(j))==1

% This is a boundary edge and a boundary integral is to be applied   

   for p=1:nip
   %sam as fortran code...
   if j==1
     eta=sqrt(3.)*0.5*(xt(p)+1);
     xi=(-1+xt(p))*(-0.5);
   elseif j==2
     eta=sqrt(3.)*(-0.5)*(xt(p)-1);
     xi=-1*(1+xt(p))*0.5;
   else
    eta=0.;
    xi=xt(p);
   end
   phx=rbasisxt(:,p+(nip*(j-1)));
   phy=rbasisyt(:,p+(nip*(j-1)));
   ph=[phx phy];
   
   J=maptri(xyt,xi,eta,localblendt);
   detJ=det(J);
   Jinv=(inv(J))';
   phxy=Jinv*ph';
   phxy=phxy';
   
   nm=inv(J)'*norml(:,j);
   nm=nm./norm(nm);
   %NB diffferent from fortran code!!!
   if j==1
     ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*0.5));
     ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*0.5));
     v=sqrt((ddgx^2)+(ddgy^2));
   elseif j==2
      ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*-0.5));
      ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*-0.5));
      v=sqrt((ddgx^2)+(ddgy^2));
   else
      ddgx=J(1,1);
      ddgy=J(2,1);
      v=sqrt((ddgx^2)+(ddgy^2));
   end   
   [xp,yp]=getcoordt(xyt,xi,eta,localblendt);
   
%   wave propagation   
%   curl=((kx^2)+(ky^2))*i*(cos((kx*xp)+(ky*yp))-i*sin((kx*xp)+(ky*yp)));
%   curl=-((kx^2)+(ky^2))*exp((kx*xp)+(ky*yp));

    %  use the problem file
    fun=probdata.neufun;
    arg=probdata.neufunarg;
    index=lbc(j);
    tangcurl=fun(xp,yp,index,nm,arg);

   for k=1:esize
     rhsel(k)=rhsel(k)-(phxy(k,1)*tangcurl(1)+phxy(k,2)*tangcurl(2))*wt(p)*v;
   end
   end
 end
end

% Now add source term for dielectric
nip=length(intxi);
if lfnum~=1
for p=1:nip
%cph=curlbasist(order,intxi(p),inteta(p));
cph=curlbast(:,p);
J=maptri(xyt,intxi(p),inteta(p),localblendt);
detJ=det(J);
cphxy=cph./detJ;

[xp,yp]=getcoordt(xyt,intxi(p),inteta(p),localblendt);


% compute curl src term
    fun=probdata.curlsrcfun;
    arg=probdata.curlsrcfunarg;
    mycurl=fun(xp,yp,arg);

   for k=1:esize
     rhsel(k)=rhsel(k)-((1./mur)-1)*cphxy(k)*mycurl*detJ*intw(p);
   end

   phx=basisxt(:,p);
   phy=basisyt(:,p);
   ph=[phx phy];
   detJ=det(J);
   if detJ<=0
   disp('Error cannot continue')
   return
   end
   Jinv=(inv(J))';
   phxy=Jinv*ph';
   phxy=phxy';

% compute src term
    fun=probdata.srcfun;
    arg=probdata.srcfunarg;
    ei=fun(xp,yp,arg);


   for k=1:esize
   for kk=1:2
     rhsel(k)=rhsel(k)+omega^2*(epr-1)*phxy(k,kk)*ei(kk)*detJ*intw(p);
   end
   end

end
end
