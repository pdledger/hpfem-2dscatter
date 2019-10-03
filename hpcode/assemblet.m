function [rhs,stiff,basisxt,basisyt,curlbast,intxi,inteta,intw,rbasisxt,rbasisyt,xt,wt,stiffabs,dirval]=assemblet(nelemt,Mesht,...
         order,map,irn,icn,nz,unkl,unkh,...
         unkti,glob,nunk,dir,nunki,edges,nelemq,x,w,rhs,dirval...
         ,xmin,ymin,tx,ty,probdata,coefft,coeffq,I,J,X,Xabs)
	 
coord=Mesht.Coordinates;
intmat=Mesht.Elements;
fnum=Mesht.Fnum;
% assemble the triangular elements performing static condensation on the fly

% compute the integration weights and locations 
[intxi,inteta,intw]=gautri(2*(order+1)+order);

% NB different set to those used with quadrilaterials! these are from -1 to +1!!!
[xt,wt]=gaussquad(-1.,1,2*(order+1)+order); 
% work out basis at each integration point
if nelemt ~= 0
[basisxt,basisyt,curlbast,rbasisxt,rbasisyt]=myevalt(...
          order,intxi,inteta,xt);
else
basisxt=[];
basisyt=[];
curlbast=[];
rbasisxt=[];
rbasisyt=[];
end


% NB omega^2 must not be an eigenvalue of Kx =lambda Mx
%ie for the box (0,pi)^2 not 1,2 etc
omega=probdata.omega;
theta=probdata.theta;
mat=probdata.mat;

[dirval]=dirichlett(edges,xt,wt,glob,coord,intmat,unkl,unkh,dir,...
                    rbasisxt,rbasisyt,nelemt,nelemq,dirval,order,probdata,coefft);

esize=3;
if order > 0
esize=(order+1)*(order+2);
end
cont=3*(order+1);

%for i=1:nelemt
for i=1:nelemt

   for j=1:3
      for k=1:2
         xyt(j,k)=coord(intmat(i,j),k);
      end
   end
   
   for j=1:3
      lbc(j)=edges(glob(i+nelemq,j),3);
   end

   for j=1:3
      for k=1:2
      localblendt(j,k)=coefft(i,j,k);
      end
   end

%  work out epr and mur in each element
   mur=mat(fnum(i),1);
   epr=mat(fnum(i),2);



% compute the elemental quadrilateral mass matrix
 [massel,nmassel]=masst(xyt,order,intxi,inteta,intw,basisxt,basisyt,xmin,ymin,tx,ty,epr,localblendt);
  
% compute the elemental quadrilateral stiffness matrix  
  [stiffel,nstiffel]=stifft(xyt,order,intxi,inteta,intw,curlbast,xmin,ymin,tx,ty,mur,localblendt);

% determine the elemental contribution to the rhs
%  rhsel=rhst(xyt,order,xt,wt,rbasisxt,rbasisyt,lbc,probdata,localblendt);
rhsel=rhst(xyt,order,xt,wt,rbasisxt,rbasisyt,lbc,probdata,localblendt,epr,mur,...
curlbast,basisxt,basisyt,intw,intxi,inteta,omega,fnum(i));

% perform static condensation
  
  [kcc,rc,kccabs]=static1(massel,stiffel,rhsel,order,omega,esize,cont,nstiffel,nmassel);
  
% begin assembly proceedure 
% create a list of unknown numbers in the order of the elemental matrix
% lowest order block
  for j=1:3
     lunk(j)=unkl(glob(i+nelemq,j));
     ldr(j)=dir(i+nelemq,j);
  end
% Higher order block
  if order> 0
     for j=1:3
     for p=1:order
        lunk(j+3*p)=unkh(glob(i+nelemq,j),p);
        ldr(j+3*p)=dir(i+nelemq,j)^(p-1);
     end
     end
%Interior block
     for j=1:order*(order-1)+order-1
        lunk(3*(order+1)+j)=unkti(i,j);
        ldr(3*(order+1)+j)=1;
     end
  end
% Now assemble matrices in sparse format
  
  help=zeros(nunki,1);
  for j=1:cont % nb changed from cont
     row=lunk(j);  
     dr=ldr(j);
     if row > 0
        for k=1:cont % nb changed from cont
           col=lunk(k);
           dc=ldr(k);
           if col > 0
             nz=nz+1;
	     len=length(X);
	     if nz > len
	       I(2*len)=0;
	       J(2*len)=0;
	       X(2*len)=0;
	       Xabs(2*len)=0;
	     end
	     I(nz)=row;
	     J(nz)=col;
	     X(nz)=kcc(j,k)*dr*dc;
	     Xabs(nz)=kccabs(j,k)*dr*dc;
            else
	      rhs(row)=rhs(row)-(kcc(j,k)*dr*dc*dirval(abs(col)));
           end
        end
	rhs(row)=rhs(row)+(dr*rc(j));
     end
  end   
end

stiff=sparse(I(1:nz),J(1:nz),X(1:nz),nunk,nunk);
stiffabs=sparse(I(1:nz),J(1:nz),Xabs(1:nz),nunk,nunk);

clear help
clear kcc
clear kccanbs
clear massel
clear stiffel
clear rc
clear rhsel
