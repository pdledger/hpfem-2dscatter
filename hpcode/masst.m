function [massel,nmassel]=masst(xy,order,intxi,inteta,intw,basisxt,basisyt,xmin,ymin,tx,ty,epr,localblendt)

nip=length(intw);
if order==0
esize=3;
else
esize=(order+1)*(order+2);
end
massel=zeros(esize);
nmassel=zeros(esize);
for p=1:nip
%ph=basist(order,intxi(p),inteta(p));

phx=basisxt(:,p);
phy=basisyt(:,p);
ph=[phx phy];
J=maptri(xy,intxi(p),inteta(p),localblendt);
detJ=det(J);
if detJ<=0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';
phxy=Jinv*ph';
phxy=phxy';

% PML TENSOR
ten=zeros(2,2);
[xp,yp]=getcoordt(xy,intxi(p),inteta(p),localblendt);
c=-log10(1e-5);
pw=5;
if abs(xp) >= xmin
z(1)=1-pw*j*c*(abs(xp)-xmin)^(pw-1)/tx^pw*xp*sign(xp)-j*c*(abs(xp)-xmin)^pw/tx^pw;
else
z(1)=1;
end
if abs(yp) >= ymin
z(2)=1-pw*j*c*(abs(yp)-ymin)^(pw-1)/ty^pw*yp*sign(yp)-j*c*(abs(yp)-ymin)^pw/ty^pw;
else
z(2)=1;
end
ten(1,1)=z(2)./z(1);	  
ten(2,2)=z(1)./z(2);

tphxy=phxy*ten;

% non vectorized
% for i=1:esize
% for jj=1:esize
% for k=1:2
% massel(i,jj)=massel(i,jj)+(intw(p)*detJ*ten(k)*phxy(i,k)*phxy(jj,k)*epr);
% nmassel(i,jj)=nmassel(i,jj)+(intw(p)*detJ*phxy(i,k)*phxy(jj,k)*epr);
% end
% end
% end

% vectorized
massel=massel+(intw(p)*detJ*tphxy*phxy'*epr);
nmassel=nmassel+(intw(p)*detJ*phxy*phxy'*epr);
 
end
end
