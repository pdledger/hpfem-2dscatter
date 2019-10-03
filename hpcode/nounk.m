function [order,unkl,unkh,unkqi,unkti,nunk,nunki,ndir,nunkpmlz,newglob,newedges,xmin,ymin,xmax,ymax,nedgepml]=nounk(ne,edges,nelemq,nelemt,glob,Mesht,Meshq,order,bctype,probdata)


% flag those elements in the PML
xmin=probdata.xmin;
xmax=probdata.xmax;
ymin=probdata.ymin;			     
ymax=probdata.ymax;
tx=xmax-xmin;
ty=ymax-ymin;
coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
intmat=Mesht.Elements;
[x,w]=gaussquad(0.,1,2*(order+1));
nip=length(x);
flagedge=zeros(ne,1);
for i=1:nelemq
   flag=0;
   for j=1:4
      for k=1:2
         xyq(j,k)=coord(intmaq(i,j),k);
      end
   end
   for p=1:nip
   for q=1:nip
   [xp,yp]=getcoord(xyq,x(p),x(q));
   if abs(xp) > xmin
   flag=1;
   elseif abs(yp) > ymin
   flag=1;
   end
   end
   end
   if flag==1
   for j=1:4
   flagedge(glob(i,j))=1;
   end
   end
end
% compute the integration weights and locations 
[intxi,inteta,intw]=gautri(2*(order+1));
nipt=length(intxi);
for i=1:nelemt
   flag=0;
   for j=1:3
      for k=1:2
         xyt(j,k)=coord(intmat(i,j),k);
      end
   end   
   localblend=-1*ones(3,2);
   for p=1:nipt
   [xp,yp]=getcoordt(xyt,intxi(p),inteta(p),localblend);
   if abs(xp) > xmin
   flag=1;
   elseif abs(yp) > ymin
   flag=1;
   end
   end
   if flag==1
   for j=1:3
   flagedge(glob(i+nelemq,j))=1;
   end
   end
end

% Now renumber the edges so that those in the PML are numbered first
newedgeno=zeros(ne,1);
newno=0;
for i=1:ne
  if flagedge(i)==1
  newno=newno+1;
  newedgeno(i)=newno;
  end
end
nedgepml=newno;
for i=1:ne
  if newedgeno(i)==0
  newno=newno+1;
  newedgeno(i)=newno;
  end
end

newedges=zeros(ne,3);
for i=1:ne
for j=1:3
newedges(newedgeno(i),j)=edges(i,j);
end
end
newglob=zeros(nelemt+nelemq,4);
for i=1:nelemq
for j=1:4
newglob(i,j)=newedgeno(glob(i,j));
end
end
for i=1:nelemt
for j=1:3
newglob(i+nelemq,j)=newedgeno(glob(i+nelemq,j));
end
end

% in this function we number the unknowns
% p=0 block

unkl=zeros(ne,1);
nunk=0;
ndir=0;
% First those in the PML (zero block)
unkl=zeros(ne,1);
for i=1:nedgepml
   if newedges(i,3)==0 | (newedges(i,3)~=0 & bctype(newedges(i,3))~=2) %newedges(i,3) ~= 2  && newedges(i,3) ~=1
      nunk=nunk+1;
      unkl(i)=nunk;
   else
      ndir=ndir+1;
      unkl(i)=-1*ndir;      
   end
end
% First those in the PML (higher block)
unkh=[];
if order > 0
unkh=zeros(ne,order);
for i=1:nedgepml
   if newedges(i,3)==0 | (newedges(i,3)~=0 & bctype(newedges(i,3))~=2) %newedges(i,3) ~= 2 && newedges(i,3) ~=1
      for j=1:order
         nunk=nunk+1;
         unkh(i,j)=nunk;	 
      end
   else
      for j=1:order
         ndir=ndir+1;
         unkh(i,j)=-1*ndir;	 
      end
   end
end
end

nunkpmlz=nunk;
disp(['There are',num2str(nunkpmlz),'unknowns in the PML']);
pause(2)
% Remaining functions
for i=nedgepml+1:ne
   if newedges(i,3)==0 | (newedges(i,3)~=0 & bctype(newedges(i,3))~=2) %newedges(i,3) ~= 2  && newedges(i,3) ~=1
      nunk=nunk+1;
      unkl(i)=nunk;
   else
      ndir=ndir+1;
      unkl(i)=-1*ndir;      
   end
end
%nunkpmlz=nunk;
disp(['There are',num2str(nunkpmlz),'unknowns  the zero and PML block']);
unkqi=[];
unkti=[];
nunki=nunk;
if order > 0
% p>0 edge block
   
for i=nedgepml+1:ne
   if newedges(i,3)==0 | (newedges(i,3)~=0 & bctype(newedges(i,3))~=2) %newedges(i,3) ~= 2 && newedges(i,3) ~=1
      for j=1:order
         nunk=nunk+1;
         unkh(i,j)=nunk;	 
      end
   else
      for j=1:order
         ndir=ndir+1;
         unkh(i,j)=-1*ndir;	 
      end
   end
end

disp(['the number of edgebased unknowns is',num2str(nunk)])

% interior block quadrilaterals
nunki=nunk;
if nelemq ~= 0
for i=1:nelemq
   for j=1:2*order*(order+1)
      nunki=nunki+1;
      unkqi(i,j)=nunki;
   end
end
end

% interior block triangles
if nelemt ~= 0
for i=1:nelemt
   for j=1:order*(order-1)+order-1
      nunki=nunki+1;
      unkti(i,j)=nunki;
   end
end
end

end

disp(['the number of unknowns are',num2str(nunki)]);
disp(['the number of knowns are',num2str(ndir)]);
