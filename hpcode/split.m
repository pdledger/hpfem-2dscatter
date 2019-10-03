function [Meshsq,Meshst,splitq,orders,splitt]=split(ne,order,Meshq,nelemq,Mesht,nelemt,dir,glob,edges,coefft,coeffq)

orders=min((order)^2,1);
coord=Meshq.Coordinates;
[npoin,i]=size(coord);
intmaq=Meshq.Elements;
intmat=Mesht.Elements;

% create number on splitedge
 k=0;
for i=1:ne
for j=1:orders+1
   k=k+1;
   splitedge(i,j)=k;
end
end
splitq=zeros(nelemq,(orders+3)^2);
      
% Quantities for triangular splittings
nintpt=0;
for i=1:orders
   nintpt=nintpt+i;
end
%  number of elements per triangle
nelpt=4;
for i=1:orders
   nelpt=nelpt+3+(2*i);
end
% number of points per triangle
nptri=6;
for i=1:orders
   nptri=nptri+(3+i);
end
intman=zeros(nelemq*((orders+2)^2),4);
intmant=zeros((nelpt*nelemt),3);            
coordn=zeros(npoin+(ne*(orders+1))+(nelemq*((orders+1)^2))...
       +(nelemt*nintpt),2);
%bsidon=zeros(nboun*(orders+2),5);
% create points on split quad...
nlayer=orders+2;      
      
epetq=splitquad(orders);
[epett,layer]=splittri(orders,nptri);

% note npoin includes number of vertices
npoinn=npoin+(ne*(orders+1));
nelemnt=0;
nelemn=0;
nbounn=0;

% loop over quads and create new elements
for i=1:nelemq

   for j=1:4
   for k=1:2
      xy(j,k)=coord(intmaq(i,j),k);
   end
   end
   
   for j=1:4
      for k=1:2
      localblend(j,k)=coeff(i,j,k);
      end
   end
      
% already have 4 conrners
   con=zeros(1,(orders+3)^2);      
   con(1)=intmaq(i,1);
   con(orders+3)=intmaq(i,2);
   con((orders+3)^2)=intmaq(i,3);
   con(((orders+3)^2)-(orders+2))=intmaq(i,4);
      
%  and also each edge
%  edge 1
  if dir(i,1) > 0
     for j=1:orders+1
        con(1+j)=splitedge(glob(i,1),j)+npoin;
     end
  else
     p=orders+1;
     for j=1:orders+1
        con(1+p+1-j)=splitedge(glob(i,1),j)+npoin;
     end      
   end
% edge 2
   if dir(i,2) > 0
      for j=1:orders+1
         con(((orders+3)^2)-j)=splitedge(glob(i,2),j)+npoin;
      end
   else
      p=orders+1;
      for j=1:orders+1
         con(((orders+3)^2)-(p+1-j))=splitedge(glob(i,2),j)+npoin;
      end
   end
% edge 3
   if dir(i,3) > 0
      for j=1:orders+1
         con((orders+3)+(j*(orders+3)))=splitedge(glob(i,3),j)+npoin;
      end
   else
      p=orders+1;
      for j=1:orders+1
         con((orders+3)+((p+1-j)*(orders+3)))=splitedge(glob(i,3),j)+npoin;
      end
   end
% edge 4
   if dir(i,4) > 0
      for j=1:orders+1       
         con(((orders+3)^2)-(orders+2)-(j*(orders+3)))=splitedge(glob(i,4),j)+npoin;
      end
   else
      p=orders+1;
      for j=1:orders+1
          con(((orders+3)^2)-(orders+2)-((p+1-j)*(orders+3)))=splitedge(glob(i,4),j)+npoin;
      end
   end
      
   for j=1:(orders+3)^2
      if con(j)==0
         npoinn=npoinn+1;
         con(j)=npoinn;
      end
   end

% store con for plotting
   for j=1:(orders+3)^2
      splitq(i,j)=con(j);
  end      
%     create connectivities
   for j=1:orders+2
   for p=1:orders+2
      nelemn=nelemn+1;
      intman(nelemn,1)=con(p+((j-1)*(orders+3)));
      intman(nelemn,2)=con(p+1+((j-1)*(orders+3)));
      intman(nelemn,3)=con(p+1+(j*(orders+3)));
      intman(nelemn,4)=con(p+(j*(orders+3)));
   end
   end     
    
   for j=1:(orders+3)^2
      [coordn(con(j),1),coordn(con(j),2)]=getcoord(xy,epetq(j,1),epetq(j,2),localblend);
   end


end            
coord=Mesht.Coordinates;
splitt=[];
% Now split triangular elements
for i=1:nelemt

   for j=1:3
   for k=1:2
      xyt(j,k)=coord(intmat(i,j),k);
   end
   end

   for j=1:3
      for k=1:2
      localblendt(j,k)=coefft(i,j,k);
      end
   end
   
% already have 3 conrners
   cont=zeros(nptri);  
      
   cont(1)=intmat(i,2);
   cont(nptri-(2+orders))=intmat(i,3);
   cont(nptri)=intmat(i,1);
      
% and also each edge
% edge 1
   if dir(i+nelemq,1)> 0
      for j=1:orders+1
         cont(layer(nlayer-j+1,1+nlayer-j))=splitedge(glob(i+nelemq,1),j)+npoin;
      end
   else
      p=orders+1;
      for j=1:orders+1
        cont(layer(j+1,j+1))=splitedge(glob(i+nelemq,1),j)+npoin;
      end      
   end
%     edge 2
   if dir(i+nelemq,2) > 0
     for j=1:orders+1
        cont(layer(j+1,1))=splitedge(glob(i+nelemq,2),j)+npoin;
     end
   else
     p=orders+1;
     for j=1:orders+1
        cont(layer(nlayer-j+1,1))=splitedge(glob(i+nelemq,2),j)+npoin;
     end
   end
%     edge 3
   if dir(i+nelemq,3) > 0
     for j=1:orders+1
        cont(layer(nlayer+1,1+j))=splitedge(glob(i+nelemq,3),j)+npoin;
     end
   else
     p=orders+1;
     for j=1:orders+1
        cont(layer(nlayer+1,1+nlayer-j))=splitedge(glob(i+nelemq,3),j)...
        +npoin;
     end
   end
      
   for j=1:nptri
      if cont(j) == 0
         npoinn=npoinn+1;
         cont(j)=npoinn;
      end
   end

%     store con data for plotting....    
   for  j=1:nptri
      splitt(i,j)=cont(j);
   end
   
   for p=1:nlayer
      for j=1:p
         nelemnt=nelemnt+1;
         intmant(nelemnt,1)=cont(layer(p-1+1,j));
         intmant(nelemnt,2)=cont(layer(p+1,j));
         intmant(nelemnt,3)=cont(layer(p+1,j+1));
      end
      if p>= 2
         for j=1:p-1
            nelemnt=nelemnt+1;
            intmant(nelemnt,1)=cont(layer(p-1+1,j+1));
            intmant(nelemnt,2)=cont(layer(p-1+1,j));
            intmant(nelemnt,3)=cont(layer(p+1,j+1));
         end
      end
   end
                
   for j=1:nptri
      [coordn(cont(j),1),coordn(cont(j),2)]=getcoordt(xyt,epett(j,1),epett(j,2),localblendt);   
   end

           
end            


if nelemn~=nelemq*((orders+2)^2)
  disp('incorrect number of elements created');
end
      
if nelemnt~=nelemt*nelpt
  disp('wrong number of triangular elements created');
end

%if nbounn~=nboun*(orders+2)
%   disp('wrong number of bounary segments created');
%end
      
if npoinn~=(ne*(orders+1))+npoin+(nelemq*((orders+1)^2))+...
   (nelemt*nintpt)
   disp('incorrect number of points created');
end      

Meshsq.Coordinates=coordn;
Meshsq.Elements=intman;

Meshst.Coordinates=coordn;
Meshst.Elements=intmant;
% call to plotting function
figure
hold on
if nelemt~=0
plot_Mesh(Meshst);
end
if nelemq~=0
plot_Mesh(Meshsq);
end
