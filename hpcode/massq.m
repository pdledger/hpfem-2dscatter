function [massel,nmassel]=massq(xy,order,x,w,basisx,basisy,xmin,ymin,tx,ty,localblend)

nip=length(x);
esize=4*(order+1)+2*order*(order+1);
massel=zeros(esize);
nmassel=zeros(esize);
for p=1:nip
for q=1:nip

phx=basisx(:,q+(nip*(p-1)));
phy=basisy(:,q+(nip*(p-1)));
ph=[phx phy];
J=mapquad(xy,x(p),x(q),localblend);
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';
phxy=Jinv*ph';
phxy=phxy';

% PML TENSOR
[xp,yp]=getcoord(xy,x(p),x(q),localblend);
if abs(xp) >= xmin
z(1)=1-j*(10*((abs(xp)-xmin)/0.5).^4*(xp*sign(xp))+...
          ((abs(xp)-xmin)/0.5)^5);
else
z(1)=1;
end
if abs(yp) >= ymin
z(2)=1-j*(10*((abs(yp)-ymin)/0.5).^4*(yp*sign(yp))+...
          ((abs(yp)-ymin)/0.5)^5);
else
z(2)=1;
end
ten(1)=z(2)./z(1);	  
ten(2)=z(1)./z(2);


for i=1:esize
for jj=1:esize
for k=1:2
massel(i,jj)=massel(i,jj)+(w(p)*w(q)*detJ*ten(k)*phxy(i,k)*phxy(jj,k));
nmassel(i,jj)=nmassel(i,jj)+(w(p)*w(q)*detJ*phxy(i,k)*phxy(jj,k));
end
end
end

end
end
