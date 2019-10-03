function [nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals]=read_Mesh(filename,rot)
%    Input
%   filename = string containing the path and the filename edge.plt
%   Output...
%   A picture of the mesh

fid=fopen(filename,'r');
n=fscanf(fid,'%5d',1)
a=fscanf(fid,'%s',11)
data=fscanf(fid,'%5d',4);
nelemq=data(1);
nelemt=data(2);
npoin=data(3);
nboun=data(4)
elementsq=[];
a=fscanf(fid,'%s',1);
% read in the quadrilateral connectivities
for i=1:nelemq
	data=fscanf(fid,'%5d',5);
	data=data';
	elementsq=[elementsq; data(2:5)];
end
% read in the triangular connectibities
elementst=[];
for i=1:nelemt
	data=fscanf(fid,'%5d',5);
	data=data';
	elementst=[elementst; data(2:4)];
end
a=fscanf(fid,'%s',1);
% read in the coordinates.
coordinates=[];
for i=1:npoin
	data=fscanf(fid,'%d %g %g %g %g',5);
	data=data';
	coordinates=[coordinates; data(2:3)];
end

% apply rotation
for i=1:npoin
x=coordinates(i,1);
y=coordinates(i,2);
coordinates(i,1)=x*cos(rot)-y*sin(rot)+1e-10;
coordinates(i,2)=x*sin(rot)+y*cos(rot)+1e-10;
end

a=fscanf(fid,'%s',1);
% read in unknowns (blank).
for i=1:npoin
	data=fscanf(fid,'%d %g %g %g %g',5);
end
a=fscanf(fid,'%s',2);
% read in boundary data
bsido=[];
for i=1:nboun
	data=fscanf(fid,'%d %d %d %d %d',5);
	data=data';
	bsido=[bsido; data(1:5)];
end
bcvals=zeros(nboun,1);
a=fscanf(fid,'%s',2);
a=fscanf(fid,'%s',3);
% read in the edge data
data=fscanf(fid,'%d %d %d',3);
ne=data(1);
kl=data(2);
% read in the element flag
a=fscanf(fid,'%s',2);
for i=1:nelemt+nelemq
	data=fscanf(fid,'%d %d',2);
end
% read in the global array
a=fscanf(fid,'%s',1);
glob=[];
for i=1:nelemt+nelemq
	data=fscanf(fid,'%d %d %d %d %d',5);
	data=data';
	glob=[glob; data(2:5)];	
end
a=fscanf(fid,'%s',1);
edges=[];
for i=1:ne
	data=fscanf(fid,'%d %d %d %d',4);
	data=data';
	edges=[edges; data(2:4)];
end
a=fscanf(fid,'%s',1);
dir=[];
for i=1:nelemq+nelemt
	data=fscanf(fid,'%d %g %g %g %g',5);
	data=data';
	dir=[dir; data(2:5)];
end


% Mesh Data
Mesht.Coordinates=coordinates;
Meshq.Coordinates=coordinates;
Mesht.Elements=elementst;
Meshq.Elements=elementsq;
fnum=ones(nelemq,1);
Meshq.Fnum=fnum;
Mesht.Fnum=[];


% call to plotting function
figure
plot_Mesh(Mesht);
plot_Mesh(Meshq);

figure
for i=1:nboun
plot(coordinates([bsido(i,1); bsido(i,2)],1),coordinates([bsido(i,1); bsido(i,2)],2),'-x');
hold on
end
hold off
