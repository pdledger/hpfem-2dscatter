function [data]=rcs(sol,kx,ky,basisx,basisy,curlbasis,nelemq,Meshq,order,unkl,unkh,...
         unkqi,unkti,glob,nunk,dir,nunki,x,w,Mesht,intxi,inteta,intw,...
	     basisxt,basisyt,curlbast,nelemt,omega,edges,rbasisx,rbasisy,...
	     probdata)
	     
% This function implements a simple version of RCS calculation. 
% There are more accurate methods for doing this, these can be
% found in the thesis
% P.D. Ledger An hp adaptive proceedure for electromagnetic scattering
% problems. Swansea 2002.

omega=probdata.omega;
theta=probdata.theta;
mat=probdata.mat;
bcscatter=probdata.bcscatter;
bcmat=probdata.bcmat;

coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
intmat=Mesht.Elements;
fnum=Meshq.Fnum;

nip=length(x);
% 4 bc edge possibilities
% Edge 1 y=0
crbasis=[];
for p=1:nip
cph=curlbasisq(order,x(p),0);
crbasis=[crbasis, cph];
end
% Edge 2 y=1
for p=1:nip
cph=curlbasisq(order,x(p),1);
crbasis=[crbasis, cph];
end
% Edge 3 x=1
for p=1:nip
cph=curlbasisq(order,1,x(p));
crbasis=[crbasis, cph];
end
% Edge 4 x=0
for p=1:nip
cph=curlbasisq(order,0,x(p));
crbasis=[crbasis, cph];
end

nonbc=0;
rcslabel=[];
for i=1:nelemq
if fnum(i)==1
for jj=1:4
if edges(glob(i,jj),3)==bcscatter
nonbc=nonbc+1;
rcslabel(nonbc,1)=i;
rcslabel(nonbc,2)=jj;
end
end
end
end

norml=[0 -1; 0 1; 1 0; -1 0]';
data=[];
dphi=pi/200.;
for npec=-200:200
phi=dphi*npec;
sigma=0;

for non=1:nonbc
   k=rcslabel(non,1);
   p=rcslabel(non,2);
% elemetal size
   esize=2*(order+1)*(order+2);

   for l=1:4
      for kk=1:2
         xyq(l,kk)=coord(intmaq(k,l),kk);
      end
   end
   
   for l=1:4
      for kk=1:2
      localblend(l,kk)=coeff(i,l,kk);
      end
   end
   
   % lowest order block
  for l=1:4
     lunk(l)=unkl(glob(k,l));
     ldr(l)=dir(k,l);
  end
% Higher order block
  if order> 0
     for l=1:4
     for pp=1:order
        lunk(l+4*pp)=unkh(glob(k,l),pp);
        ldr(l+4*pp)=dir(k,l)^(pp-1);
     end
     end
%Interior block
     for l=1:2*order*(order+1)
        lunk(4*(order+1)+l)=unkqi(k,l);
        ldr(4*(order+1)+l)=1;
     end
   end
  
   int=0;
   for kk=1:nip
      
      if p==1
      xi=x(kk);     
      eta=0;
      elseif p==2
      xi=x(kk);
      eta=1;
      elseif p==3
      xi=1;
      eta=x(kk);
      else
      xi=0;
      eta=x(kk);
      end

      J=mapquad(xyq,xi,eta,localblend);
      detJ=det(J);
      Jinv=(inv(J))';
      % Basis
%      ph=basisq(order,xi,eta);
      phx=rbasisx(:,kk+(nip*(p-1)));
      phy=rbasisy(:,kk+(nip*(p-1)));
      ph=[phx phy];
      phxy=Jinv*ph';
      phxy=phxy';
      %Curl Basis
%      cph=curlbasisq(order,xi,eta);
      cph=crbasis(:,kk+(nip*(p-1)));
      cphxy=cph./detJ;
      
      etil=[0;0];
      cetil=0;
      for kkk=1:esize
      if lunk(kkk) ~=0
        etil(1)=etil(1)+(phxy(kkk,1).*sol(lunk(kkk)).*ldr(kkk));      
        etil(2)=etil(2)+(phxy(kkk,2).*sol(lunk(kkk)).*ldr(kkk));
	    cetil=cetil+(cphxy(kkk).*sol(lunk(kkk)).*ldr(kkk));
      end
      end 
      mag=-1.*cetil./(j*omega);

      nm=inv(J)'*norml(:,p);
      nm=nm./norm(nm);
   
      if p==1 || p==2
      v=sqrt(J(1,1)^2+J(2,1)^2);
      else
      v=sqrt(J(1,2)^2+J(2,2)^2);
      end

      val(3)=(nm(2)*etil(1))-(nm(1)*etil(2));
      val(2)=-nm(1)*mag;
      val(1)=nm(2)*mag;

     [xp,yp]=getcoord(xyq,xi,eta,localblend);
      solc=-(val(1)*(sin(phi)))+(val(2)*(cos(phi)))+val(3);
      fact2=(omega)*((xp*(cos(phi)))+(yp*(sin(phi))));
      fact3=cos(fact2)+j*sin(fact2);

      int=int+(fact3*v*solc*w(kk));

      end

   sigma=sigma+int;
end
val=(omega/4)*((abs(sigma))^2);
val=10.*(log10(val));
data=[data;(phi*180.)/pi,val];
end
figure
plot(data(:,1),data(:,2))
