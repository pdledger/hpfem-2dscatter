
function probdata=problem1()
% PEC NACA

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to H(curl) elements with const. Tang. Component.
order=0;
probdata.order=order;

% define basic parameters
tetm=2;           % TE =1 TM =2
% note that the wavenumber kappa is called omega in the code
omega=2*pi;      % =K_0 since epr=1, mur=1
theta=0;      % angle of incidance in radians
pmlt=0.5;        % PML Thickness note that the domain size must be large 
                 % enough to accomodate this thickness
probdata.tetm=tetm;
probdata.omega=omega;
probdata.theta=theta;

% define the materials
mu0=1;
ep0=1;
mat0=[mu0 ep0];
probdata.mat=mat0;

% curved geometry information
% lingeom = 0 Curved
% lingeom = 1 Linear
lingeom=1; 
rin=1;
probdata.rin=rin;
probdata.lingeom=lingeom;

% set mesh spacing
h=0.2;
probdata.h=h;

% Define problem geometry
% Define problem geometry
node=[1 1; 
      0 1;
      0 0;
      1 0];
% centre airfoil
node(:,1)=node(:,1)-0.5;
node(:,2)=node(:,2)-0.5;
[np dum]=size(node);
cnect=[];
for i=1:np-1
cnect  = [ cnect; i i+1];
end
cnect=[cnect;np 1]
[nsegscat,dum]=size(cnect);
scatterer=node


node   = [node; 
         -2 -2; 
          2 -2;
          2  2;
         -2 2];  
cnectpml=[ np+1 np+2; np+2 np+3; np+3 np+4; np+4 np+1];   
[nsegpml dum]=size(cnectpml);
cnect  = [cnect; cnectpml]
probdata.node=node;
probdata.cnect=cnect;
probdata.face=[];


% flag those elements in the PML
xmin=2-pmlt;
xmax=2;
ymin=2-pmlt;			     
ymax=2;
probdata.xmin=xmin;
probdata.xmax=xmax;
probdata.ymin=ymin;
probdata.ymax=ymax;


% use mesh2d and not meshfaces
probdata.meshfaces=0;

% persibe function for mesh refinement
arg={h};
probdata.meshfun=@meshref
probdata.meshfunarg=arg;


% define boundary data
% bctype
% bctype 1 Neumann type
% bctype 2 Dirichlet

% bcflag
% The flag to denote the individual boundary segments
% each boundary segment must have a flag and a type!
%        scatterer          pml
bcflags=[1*ones(nsegscat,1);2*ones(nsegpml,1)];
if tetm==1
bctype= [2;  2];
else
bctype= [1;  2];
end
% ie each boundary is Dirichlet type and there are different Dirichlet functions
probdata.bcflags=bcflags;
probdata.bctype=bctype;

% define the function handle for Dirichlet BC's
arg.theta=theta;
arg.omega=omega;
probdata.dirfun=@problem1dir
probdata.dirfunarg=arg;

% define the function handle for Neumann BC's
probdata.neufun=@problem1neu
probdata.neufunarg=arg;

% define the function handle for the source terms
probdata.curlsrcfun=@problem1curlsrc
probdata.curlsrcfunarg=arg;

probdata.srcfun=@problem1src
probdata.srcfunarg=arg;

% define the contour and material for computing the RCS
bcscatter=1;
bcmat=1;
exactrcs=0; %exactrcs=1 available otherwise 0
probdata.bcscatter=bcscatter;
probdata.bcmat=bcmat;
probdata.exactrcs=exactrcs;
if exactrcs==1
arg.omega=omega;
arg.a=rin;
arg.tetm=tetm;
arg.theta=theta;
probdata.rcsfun=@rcsexact
probdata.rcsfunarg=arg
end


function bcval=problem1dir(x,y,index,nm,arg)
theta=arg.theta;
omega=arg.omega;

if index==1
% scatterer E^i
   e=[-sin(theta); cos(theta)]*(cos(omega*((cos(theta)*x)+(sin(theta)*y)))-...
   j*sin(omega*((cos(theta)*x)+(sin(theta)*y))));;

% n x E^s = - n x E^ i
bcval = - (nm(1)*e(2)-nm(2)*e(1));     

elseif index==2
% far field 
bcval=0;
else
disp('Not defined for this index')
end


function bcval=problem1neu(x,y,index,nm,arg)
theta=arg.theta;
omega=arg.omega;

if index==1
% scatterer curl E^i
     curle=-i*omega*(cos(omega*((x*cos(theta))...
     +(y*sin(theta))))+...
     -i*sin(omega*((x*cos(theta))+(y*sin((theta))))));

% n x curl E^s = - n x curl E^ i
bcval = - [nm(2)*curle; -nm(1)*curle];     
     
else
disp('Not defined for this index')
end



function src=problem1curlsrc(x,y,arg)
theta=arg.theta;
omega=arg.omega;

src=-i*omega*(cos(omega*((x*cos(theta))...
     +(y*sin(theta))))+...
     -i*sin(omega*((x*cos(theta))+(y*sin((theta))))));

function src=problem1src(x,y,arg)
theta=arg.theta;
omega=arg.omega;

src=[-sin(theta); cos(theta)]*(cos(omega*((cos(theta)*x)+(sin(theta)*y)))-...
   j*sin(omega*((cos(theta)*x)+(sin(theta)*y))));;


function rcsval=rcsexact(phi,arg)
tetm=arg.tetm;
omega=arg.omega;
a=arg.a;
theta=arg.theta;

if tetm==1
%exact tercs
rcsval=teexactrcs(omega,a,phi,theta);
else
%exact tmrcs
rcsval=tmexactrcs(omega,a,phi,theta);
end

function h = meshref(x,y,arg)
hglob=arg;

% User defined size function for uniform spacing modify to allow non-uniform
% spacing
r=sqrt(x.^2+y.^2);
h=hglob*sqrt(10*r);
[m n]=size(r);
minh=1e-4*hglob*ones(m,n);
h=max(h,minh);

h=hglob*ones(size(x));
