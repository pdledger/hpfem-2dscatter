function [sol]=assembleqint(nelemq,Meshq,order,unkl,unkh,unkqi,...
        unkti,glob,nunk,dir,nunki,edges,x,w,sol,basisx,basisy,curlbasis,...
        rbasisx,rbasisy,xmin,ymin,tx,ty,probdata,coefft,coeffq)

coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
% assemble the quadrilateral elements performing static condensation on the fly

% NB omega^2 must not be an eigenvalue of Kx =lambda Mx
%ie for the box (0,pi)^2 not 1,2 etc
omega=probdata.omega;
theta=probdata.theta;

esize=4*(order+1)+2*order*(order+1);
cont=4*(order+1);
inter=2*order*(order+1);

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
  
  for j=1:inter
     row=lunk(j+cont);  
     dr=ldr(j+cont);
     if row ~= 0
% create a list of which columns relate to which nonz enteries
        sol(row)=lvali(j)*dr;
     end
  end   
end




