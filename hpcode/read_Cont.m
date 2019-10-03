function read_Cont(filename)
%    Input
%   filename = string containing the path and the filename cont.plt 
%   Output...
%   A contour plot of the vairables

fid=fopen(filename,'r');
n=fscanf(fid,'%5d',1);
%a=fscanf(fid,'%s',11);
data=fscanf(fid,'%5d',4);
nelemq=data(1);
nelemt=data(2);
npoin=data(3);
elementsq=[];
%a=fscanf(fid,'%s',1);
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
%a=fscanf(fid,'%s',1);
% read in the coordinates.
coordinates=[];
for i=1:npoin
	data=fscanf(fid,'%d %g %g %g %g',5);
	data=data';
	coordinates=[coordinates; data(2:3)];
end
% read in the values
values=[];
for i=1:npoin
	data=fscanf(fid,'%d %g %g %g %g',5);
	data=data';
	values=[values; data(2:3)];
end


Mesht.Coordinates=coordinates;
Meshq.Coordinates=coordinates;
Mesht.Elements=elementst;
Meshq.Elements=elementsq;
% call to plotting function
figure
plot_Mesh(Mesht);
plot_Mesh(Meshq);
% fix number of output values to 2
nout=2;

for i=1:nout
figure
% call to plotting of values
plot_LFE(values(:,i),Mesht)
hold on;
plot_LFE(values(:,i),Meshq)
colorbar
hold off;
end

