function [coefft,coeff,Mesht,Meshq]=getcoeff(nelemq,nelemt,Mesht,Meshq,edges,glob,ne,coefft,coeff,bcscatter,rin);

coord=Mesht.Coordinates;
intmat=Mesht.Elements;
intmaq=Meshq.Elements;

for i=1:nelemq
for j=1:4
if edges(glob(i,j),3)==bcscatter
% this is a scatterer edge
% set values of coeff
n1=edges(glob(i,j),1);
n2=edges(glob(i,j),2);
coeff(i,j,1)=sqrt(coord(n1,1)^2+coord(n1,2)^2);
coeff(i,j,2)=sqrt(coord(n2,1)^2+coord(n2,2)^2);
end
end
end

for i=1:nelemt
for j=1:3
if edges(glob(i+nelemq,j),3)==bcscatter
% this is a scatterer edge
% set values of coeff
if j~=3
n1=intmat(i,j);
n2=intmat(i,j+1);
else
n1=intmat(i,j);
n2=intmat(i,1);
end
r1=sqrt(coord(n1,1)^2+coord(n1,2)^2);
coord(n1,1)=rin*coord(n1,1)/r1;
coord(n1,2)=rin*coord(n1,2)/r1;
r2=sqrt(coord(n2,1)^2+coord(n2,2)^2);
coord(n2,1)=rin*coord(n2,1)/r2;
coord(n2,2)=rin*coord(n2,2)/r2;

coefft(i,j,1)=sqrt(coord(n1,1)^2+coord(n1,2)^2);
coefft(i,j,2)=sqrt(coord(n2,1)^2+coord(n2,2)^2);
end
end
end
Mesht.Coordinates=coord;
Meshq.Coordinates=coord;
