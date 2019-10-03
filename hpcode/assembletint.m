function [sol]=assembletint(nelemt,Mesht,...
         order,unkl,unkh,unkti,glob,nunk,dir,nunki,edges,nelemq,xt,wt,...
	 basisxt,basisyt,curlbast,intxi,inteta,intw,rbasisxt,rbasisyt,sol,xmin,ymin...
         ,tx,ty,probdata,coefft,coeffq)

coord=Mesht.Coordinates;
intmat=Mesht.Elements;
fnum=Mesht.Fnum;

% assemble the triangular elements performing static condensation on the fly

% NB omega^2 must not be an eigenvalue of Kx =lambda Mx
%ie for the box (0,pi)^2 not 1,2 etc
omega=probdata.omega;
theta=probdata.theta;
mat=probdata.mat;

esize=3;
inter=0;
if order > 0
esize=(order+1)*(order+2);
inter=order*(order-1)+order-1;
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

% create a (local) list of found values
  lval=zeros(cont,1);
  for j=1:cont
  if lunk(j) ~= 0
  lval(j)=sol(lunk(j))*ldr(j);
  else
  lval(j)=0;
  end
  end

% perform static condensation
  [lvali]=static2(massel,stiffel,rhsel,order,omega,lval,esize,cont);
  
% begin assembly proceedure 
  if inter~=0
  for j=1:inter
     row=lunk(j+cont);  
     dr=ldr(j+cont);
     if row ~= 0
% create a list of which columns relate to which nonz enteries
        sol(row)=lvali(j)*dr;
     end
  end
  end   
end

