
%create inner circle of points
np=20;
dphi=pi/np;
theta  = (-pi:dphi:(pi-dphi))';
node   = [cos(theta) sin(theta)];
%include outer surface
node   = [node; -2 -2; 
                 2 -2;
		 2  2;
		 -2 2];
cnect  = [ 1:(2*np); 2:(2*np) 1]';
cnect  = [cnect; 2*np+1 2*np+2; 2*np+2 2*np+3; 2*np+3 2*np+4; 2*np+4 2*np+1];		 
%create mesh
[p,t]=mesh2d(node,cnect);

%creare the list of boundaries
[nelem dum]=size(t);
[npoin dum]=size(p);

[bsido]=my_data(p,t)
