function [rhs,basisx,basisy,curlbasis,x,w,rbasisx,rbasisy,dirval,tx,ty,I,J,X,Xabs,nz]=assembleq(nelemq,Meshq,...
         order,map,irn,icn,nz,unkl,unkh,unkqi,...
         unkti,glob,nunk,dir,nunki,edges,xmin,ymin,xmax,ymax,probdata,coefft,coeff)

coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
% assemble the quadrilateral elements performing static condensation on the fly

% compute the integration weights and locations
[x,w]=gaussquad(0.,1,2*(order+1));

% work out basis at each integration point
if nelemq~=0
[basisx,basisy,curlbasis,rbasisx,rbasisy]=myeval(order,x,w);
else
basisx=[];
basisy=[];
curlbasis=[];
rbasisx=[];
rbasisy=[];
end

% zero mass and stiffness matrices for assembly
I=zeros(nunk,1);
J=zeros(nunk,1);
X=zeros(nunk,1);
Xabs=zeros(nunk,1);
rhs=zeros(nunk,1);
nz=0;

% NB omega^2 must not be an eigenvalue of Kx =lambda Mx
%ie for the box (0,pi)^2 not 1,2 etc
omega=probdata.omega;
theta=probdata.theta;

% Apply Dirichlet Boundary Conditions
[dirval]=dirichletq(edges,x,w,glob,coord,intmaq,unkl,unkh,dir,...
                             rbasisx,rbasisy,nelemq,order,coeff);

tx=xmax-xmin;
ty=ymax-ymin;
% elemental matrix size
esize=4*(order+1)+2*order*(order+1);
cont=4*(order+1);

for i=1:nelemq

   for j=1:4
      for k=1:2
         xyq(j,k)=coord(intmaq(i,j),k);
      end
   end
   
   for j=1:4
      lbc(j)=edges(glob(i,j),3);
   end

   for j=1:4
      for k=1:2
      localblend(j,k)=coeff(i,j,k);
      end
   end


% compute the elemental quadrilateral mass matrix
  [massel,nmassel]=massq(xyq,order,x,w,basisx,basisy,xmin,ymin,tx,ty,localblend);
  
  
% compute the elemental quadrilateral stiffness matrix  
  [stiffel,nstiffel]=stiffq(xyq,order,x,w,curlbasis,xmin,ymin,tx,ty,localblend);
  

% determine the elemental contribution to the rhs
  rhsel=rhsq(xyq,order,x,w,rbasisx,rbasisy,lbc,omega,probdata,localblend);

% perform static condensation
 [kcc,rc,kccabs]=static1(massel,stiffel,rhsel,order,omega,esize,cont,nstiffel,nmassel);
  
% begin assembly proceedure 
% create a list of unknown numbers in the order of the elemental matrix
% lowest order block
  for j=1:4
     lunk(j)=unkl(glob(i,j));
     ldr(j)=dir(i,j);
  end
% Higher order block
  if order> 0
     for j=1:4
     for p=1:order
        lunk(j+4*p)=unkh(glob(i,j),p);
        ldr(j+4*p)=dir(i,j)^(p-1);
     end
     end
%Interior block
     for j=1:2*order*(order+1)
        lunk(4*(order+1)+j)=unkqi(i,j);
        ldr(4*(order+1)+j)=1;
     end
  end
% Now assemble matrices in sparse format
  
  help=zeros(nunki,1);
  for j=1:cont
     row=lunk(j);  
     dr=ldr(j);
     if row > 0
        for k=1:cont
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

clear help
clear kcc
clear kccanbs
clear massel
clear stiffel
clear rc
clear rhsel


