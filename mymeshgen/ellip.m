function [Meshq,iboun,nbo,Mesht,nelemt,nelemq,bcvals]=ellip(h,probdata)

% bcval (1) =  Neumann boundary condition value on edges of type 1
% bcval (2) =  Dirichlet boundary condition value on edges of type 2
% bcval (3) =  Dirichlet boundary condition value on edges of type 3
% bcflag = set the boundary condition type on each boundary edge

node=probdata.node;
cnect=probdata.cnect;
bcflag=probdata.bcflags;
bctype=probdata.bctype;
face=probdata.face;
bcvals=[];

for k=1:4
x(k)=node(k,1);
y(k)=node(k,2);
end

% suggest that the following is not modified.
imax=double(int16(1./h))
jmax=imax


ZERO = 0;
VXNO = 1;
rmax = 10;
TOLER  = 0.000001;
ITER = 0;

for i = 1: imax
for j = 1: jmax
	    x1(i,j) = 0.0;
	    x2(i,j) = 0.0;
	    y1(i,j) = 0.0;
	    y2(i,j) = 0.0;
end
end

% Create boundary points

spcex = x(2)-x(1);
spcey = y(2)-y(1);
deltax = spcex/real(imax-1);
deltay = spcey/real(imax-1);
j=1;
for i = 1: imax
x1(i,j) = x(1)+deltax*(i-1);
y1(i,j) = y(1)+deltay*(i-1);
x2(i,j) = x1(i,j);
y2(i,j) = y1(i,j);
end
delxsi  = deltax;
k = 0;
spcex = x(3)-x(4);
spcey = y(3)-y(4);
deltax = spcex/(imax-1);
deltay = spcey/(imax-1);

j = jmax;
for i= 1:imax
x1(i,j) = x(4)+deltax*k;
y1(i,j) = y(4)+deltay*k;
x2(i,j) = x1(i,j);
y2(i,j) = y1(i,j);
k = k +1;
end

spcex = x(3)-x(2);
spcey = y(3)-y(2);
deltax = spcex / real(jmax-1);
deltay = spcey / real(jmax-1);
delnet = deltay;
k = 0;
i= imax;
for j = 1:jmax
x1(i,j) = x(2)+deltax*k;
y1(i,j) = y(2)+deltay*k;
x2(i,j) = x1(i,j);
y2(i,j) = y1(i,j);
k = k +1;
end
spcex = x(4)-x(1);
spcey = y(4)-y(1);
deltax = spcex /real(jmax-1);
deltay = spcey / real(jmax-1);
k = 0;
i = 1;
for j = 1: jmax
x1(i,j) = x(1)+deltax*k;
y1(i,j) = y(1)+deltay*k;
x2(i,j) = x1(i,j);
y2(i,j) = y1(i,j);
k = k +1;
end

while rmax > TOLER & ITER<10000

for i = 2:imax-1

for j = 2:jmax-1

x2(i,j) = (0.25)*(x1(i+1,j)+x1(i-1,j)+x1(i,j+1)+x1(i,j-1));
y2(i,j) = (0.25)*(y1(i+1,j)+y1(i-1,j)+y1(i,j+1)+y1(i,j-1));

end

end

rmax =abs(x1(1,1)-x2(1,1));
for i = 1:imax
for j = 1:jmax
rtemp = abs(x1(i,j)-x2(i,j));
if rtemp> rmax
   rmax = rtemp;
end
end
end
for i = 1:imax
for j = 1:jmax
rtemp = abs(y1(i,j)-y2(i,j));
if rtemp> rmax
rmax=rtemp;
end
end
end
for j = 2: jmax-1
for i = 2: imax-1
x1(i,j) = x2(i,j);
y1(i,j) = y2(i,j);
end
end

ITER =ITER+1;
end

disp(['Converged in',num2str(ITER),' iterations']);

for i= 1: imax-1
  iboun(1,i) = i;
  iboun(2,i) = i+1;
  iboun(3,i)= i;
  iboun(4,i)=bcflag(1);
  iboun(5,i)=bcflag(1);
end
	nbo = imax-1;

for i = 1: imax-1
   kpt = (imax)*(jmax-1);
   nbo = nbo + 1;
   iboun(2,nbo) = kpt+(i);
   iboun(1,nbo) = iboun(2,nbo)+ 1;
   iboun(3,nbo) = (imax-1)*(jmax-1)-(imax-1)+i;
  iboun(4,nbo)=bcflag(3);
  iboun(5,nbo)=bcflag(3);
end

for j = 1: jmax-1
   nbo = nbo + 1;
   iboun(2,nbo) = ((j-1)*imax)+1;
   iboun(1,nbo) = iboun(2,nbo)+ imax;
   iboun(3,nbo) = (j-1)*(imax-1)+1;
  iboun(4,nbo)=bcflag(2);
  iboun(5,nbo)=bcflag(2);
end
for j = 1: jmax-1
   nbo = nbo + 1;
   iboun(1,nbo) = (j-1)*imax+imax;
   iboun(2,nbo) = iboun(1,nbo)+imax;
   iboun(3,nbo) = (j-1)*(imax-1) + imax-1;
  iboun(4,nbo)=bcflag(4);
  iboun(5,nbo)=bcflag(4);
end

	

ELNO = (imax-1)*(jmax-1);
ELNO  = 1;
for J = 1: jmax-1
for I = 1: imax-1
	ELCON((ELNO-1)*4+1) = I + (J-1)*imax ;
	ELCON((ELNO-1)*4+2) = I+1+(J-1)*imax;
	ELCON((ELNO-1)*4+3) = I+1+J*imax;
	ELCON((ELNO-1)*4+4) = I + J*imax;
	ELNO = ELNO + 1;
end
end

	ELNO = ELNO -1;


nvx = jmax*imax;
npoin=nvx;
nelemq=ELNO;
intmaq=[];
for I = 1:4:(ELNO)*4
intmaq=[intmaq; ELCON(I) ELCON(I+1) ELCON(I+2) ELCON(I+3) ];
end

for j = 1: jmax
for i = 1: imax
coord((i+(j-1)*imax),1)=x1(i,j);
coord((i+(j-1)*imax),2)=y1(i,j);
end
end


Mesht.Coordinates=coord;
Meshq.Coordinates=coord;
Mesht.Elements=[];
Meshq.Elements=intmaq;
fnum=ones(nelemq,1);
Meshq.Fnum=fnum;
Mesht.Fnum=[];


iboun=iboun';
nelemt=0;
figure
plot_Mesh(Mesht);
plot_Mesh(Meshq);


