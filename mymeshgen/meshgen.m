function [nelem,nelemq,Mesht,Meshq,bsido,nboun,bcvals,bctype]=meshgen(Mesht,hglob,probdata)

% An interface to the 2D unstructured meshing program written by Darren Engwirda
% outputs: 
% nelem: Number of triangles
% nelemq: Number of quadrilaterials (0)
% Mesht: Structure containing coordinates and connectivities for triangular mesh
% Meshq: as above for quadrilaterials
% bsido: Boundary condition data
% nboun: Number of boundary conditions
% bcvals: Boundary condition values

node=probdata.node;
cnect=probdata.cnect;
bcflags=probdata.bcflags;
bctype=probdata.bctype;
face=probdata.face;

bcvals =[ 0; 0; 1];

% decide which mesh generator to use
mesh2dmeshfaces=probdata.meshfaces;
hdata = [];
%hdata.fun = @h1;
%hdata.args = {hglob};
%options.dhmax = hglob;

hdata.fun = probdata.meshfun;
hdata.args = probdata.meshfunarg;
options.dhmax = hglob;


if mesh2dmeshfaces==0
[p,t]=mesh2d(node,cnect,hdata,options);
[nelem dum]=size(t);
fnum=ones(nelem,1);
else
[p,t,fnum,stats] = meshfaces(node,cnect,face,hdata,options);
end

%create the list of boundaries
% do not modifiy this code
%[bsido,p]=my_data(p,t,cnect,bcflags,node);
[bsido,p]=my_data_axi(p,t,cnect,bcflags,node);
[nboun dum]=size(bsido);

Mesht.Coordinates=p;
Mesht.Elements=t;
Mesht.Fnum=fnum;
Meshq.Coordinates=p;
Meshq.Elements=[];
Meshq.Fnum=[];
nelemq=0;
[nelem dum]=size(t);
plot_Mesh(Mesht);

