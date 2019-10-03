function probdata=problem2()
% dielectric cylinder

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to H(curl) elements with const. Tang. Component.
order=3;
probdata.order=order;

% define basic parameters
tetm=2;           % TE =1 TM =2
% note that the wavenumber kappa is called omega in the code
omega=2*pi;       % =K_0 since epr=1, mur=1
theta=pi/2;          % angle of incidance in radians
pmlt=0.5;        % PML Thickness
probdata.tetm=tetm;
probdata.omega=omega;
probdata.theta=theta;

% define the materials
mu0=1;
ep0=1;
ep1=1;
mu1=2.56
if tetm==1
mat=[mu0 ep0; mu1 ep1];
else
mat=[ep0 mu0; ep1 mu1];
end
probdata.mat=mat;


% curved geometry information
% lingeom = 0 Curved
% lingeom = 1 Linear
lingeom=0; 
rin=1;
probdata.rin=rin;
probdata.lingeom=lingeom;

% set mesh spacing
h=0.2;
probdata.h=h;

% Define the GF factor (ratio of spacing in PML to that is FS)
GF=2;

% Define problem geometry
np=20;
dphi=pi/np;
mytheta  = (pi-dphi:-dphi:-pi)';
node   = [rin*cos(mytheta) rin*sin(mytheta)];
cnect  = [ 1:(2*np); 2:(2*np) 1]';
[nsegscat dum]=size(cnect);
scatterer=node;

%include outer surface
domsize=2;
node   = [node; 
         -domsize -domsize; 
          domsize -domsize;
          domsize  domsize;
         -domsize domsize];  
cnectpml=[ 2*np+1 2*np+2; 2*np+2 2*np+3; 2*np+3 2*np+4; 2*np+4 2*np+1];   
[nsegpml dum]=size(cnectpml);
cnect  = [cnect; cnectpml];

% subdomain information
face{1}=[1:nsegpml+nsegscat];
face{2}=[1:nsegscat];

probdata.node=node;
probdata.cnect=cnect;
probdata.face=face;

% flag those elements in the PML
xmin=domsize-pmlt;
xmax=domsize;
ymin=domsize-pmlt;			     
ymax=domsize;
probdata.xmin=xmin;
probdata.xmax=xmax;
probdata.ymin=ymin;
probdata.ymax=ymax;

% use meshfaces and not mesh2d
probdata.meshfaces=1;

% persibe function for mesh refinement
arg={h,GF,domsize,pmlt};
probdata.meshfun=@meshref
probdata.meshfunarg=arg;


% define boundary data

% bctype
% bctype 1 Neumann type
% bctype 2 Dirichlet
% bctype 3 Interface (no BC's)

% bcflag
% The flag to denote the individual boundary segments
% each boundary segment must have a flag and a type!
%        scatterer          pml
bcflags=[1*ones(nsegscat,1);2*ones(nsegpml,1)];

if tetm==1
bctype= [3;  2];
else
bctype= [3;  2];
end

% For every subdomain specify the number of interfaces
% subdomain 1 region between cylinder and box
% sub domin 2 inner cylinder
nint=[1;1];
bcint=[1; 1];
probdata.nint=nint;
probdata.bcint=bcint;


probdata.bcflags=bcflags;
probdata.bctype=bctype;

% define the function handle for Dirichlet BC's
clear arg
arg.theta=theta;
arg.omega=omega;
arg.tetm=tetm;
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
exactrcs=1; %exactrcs=1 available otherwise 0
probdata.bcscatter=bcscatter;
probdata.bcmat=bcmat;
probdata.exactrcs=exactrcs;
if exactrcs==1
arg.omega=omega;
arg.a=rin;
arg.tetm=tetm;;
arg.ep1=ep1;
arg.mu1=mu1;
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
tetm=arg.tetm;

src=-i*omega*(cos(omega*((x*cos(theta))...
     +(y*sin(theta))))+...
     -i*sin(omega*((x*cos(theta))+(y*sin((theta))))));

function src=problem1src(x,y,arg)
theta=arg.theta;
omega=arg.omega;
tetm=arg.tetm;

src=[-sin(theta); cos(theta)]*(cos(omega*((cos(theta)*x)+(sin(theta)*y)))-...
   j*sin(omega*((cos(theta)*x)+(sin(theta)*y))));;


function rcsval=rcsexact(phi,arg)
tetm=arg.tetm;
omega=arg.omega;
a=arg.a;
ep1=arg.ep1;
mu1=arg.mu1;
theta=arg.theta;

if tetm==1
%exact tercs
rcsval=tedexactrcs(omega,a,phi,ep1,mu1,theta);
else
%exact tmrcs
rcsval=tmdexactrcs(omega,a,phi,ep1,mu1,theta);
end




function h = meshref(x,y,hglob,GF,domsize,pmlt)
np=length(x);
hpml=hglob/GF;

h=ones(size(x));
for i=1:np
xp=x(i);
yp=y(i);

if xp >= domsize-pmlt | xp <= -domsize+pmlt
h(i)=hpml;
elseif yp >= domsize-pmlt | yp <= -domsize+pmlt
h(i)=hpml;
else
h(i)=hglob;
end

end
