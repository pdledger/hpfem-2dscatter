function [data,res,iter]=main(h,order)
%Edge element solution of 2D scattering problems

close all

% add this directory and all subfolders
addpath(genpath('./'))

% select problem type
problem=0;

% define all data for this problem
probdata=problem3();

% Set order of elements;
order=probdata.order;

% set mesh spacing
h=probdata.h;


% generate mesh using meshge
Mesht=[];
[nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals,bctype]=meshgen(Mesht,h,problem,probdata);



% choose either 2 for TE or 3 for TM
bcscatter=2
omega=2*pi;

% Create mesh using mesh2d
Mesht=[];
[nelemt,nelemq,Mesht,Meshq,bsido,nboun]=meshgen(Mesht,bcscatter,h);


% Determine the number of edges
[ne,edges,glob,dir]=edgeno(Mesht,Meshq,bsido,nelemt,nelemq,nboun);

% Number the unknowns
[order,unkl,unkh,unkqi,unkti,nunk,nunki,ndir,nunkpmlz,glob,edges,xmin,ymin,xmax,ymax,nedgepml]=nounk(ne,edges,nelemq,nelemt,glob,Mesht,Meshq,order);

irn=[];
icn=[];
nz=0;
map=[];

% Assemble Quadrilateral elements
disp('begin assembly');
[rhs,gstiff,basisx,basisy,curlbasis,kx,ky,x,w,rbasisx,rbasisy,gstiffabs,dirval,omega,theta,xmin,ymin,tx,ty]=assembleq(nelemq,Meshq,...
         order,map,irn,icn,nz,unkl,unkh,unkqi,...
         unkti,glob,nunk,dir,nunki,edges,xmin,ymin,xmax,ymax,omega);
% Assemble Triangular Elements
disp('begin assembly');
[rhs,stiff,basisxt,basisyt,curlbast,intxi,inteta,intw,rbasisxt,rbasisyt,xt,wt,stiffabs,dirval]=assemblet(nelemt,Mesht,...
         order,map,irn,icn,nz,unkl,unkh,...        
         unkti,glob,nunk,dir,nunki,edges,kx,ky,nelemq,x,w,gstiff,rhs,gstiffabs,dirval,omega,theta,xmin,ymin,tx,ty);	 
% Solve for Edge based Degree's of freedom
disp('Solve system');
clear irn;
clear icn;
clear map;
clear gstiff;
clear gstiffabs;
sol=stiff\rhs;
res=[];

% introduce Dirichlet values back to solution vector
[sol,unkl,unkh,nunki]=introd(sol,dirval,unkl,unkh,order,ne,nunki);

disp('Assemble interiors (quads)');
 [sol]=assembleqint(nelemq,Meshq,order,unkl,unkh,unkqi,...
        unkti,glob,nunk,dir,nunki,edges,x,w,sol,kx,ky,basisx,basisy,curlbasis,...
        rbasisx,rbasisy,omega,theta,xmin,ymin,tx,ty);
disp('Assemble interiors (triangles)');
if order >1
[sol]=assembletint(nelemt,Mesht,...
         order,unkl,unkh,unkti,glob,nunk,dir,nunki,edges,kx,ky,nelemq,xt,wt,...
	 basisxt,basisyt,curlbast,intxi,inteta,intw,rbasisxt,rbasisyt,sol,...
	 omega,theta,xmin,ymin,tx,ty);
end

disp('Determine RCS');
if nelemt==0
[data]=rcs(sol,kx,ky,basisx,basisy,curlbasis,nelemq,Meshq,order,unkl,unkh,...
         unkqi,unkti,glob,nunk,dir,nunki,x,w,Mesht,intxi,inteta,intw,...
	     basisxt,basisyt,curlbast,nelemt,omega,edges,rbasisx,rbasisy,...
	     bcscatter);
else
[data]=rcst(sol,kx,ky,basisx,basisy,curlbasis,nelemq,Meshq,order,unkl,unkh,...
         unkqi,unkti,glob,nunk,dir,nunki,xt,wt,Mesht,intxi,inteta,intw,...
	     basisxt,basisyt,curlbast,nelemt,omega,edges,rbasisxt,rbasisyt,...
	     bcscatter)	     
end
[Meshsq,Meshst,splitq,orders,splitt]=split(ne,order,Meshq,nelemq,Mesht,nelemt,dir,glob,edges);

mycontour(nelemt,nelemq,orders,Meshst,Meshsq,order,splitq,unkl,unkh...
         ,unkqi,unkti,dir,glob,Meshq,Mesht,sol,splitt);
