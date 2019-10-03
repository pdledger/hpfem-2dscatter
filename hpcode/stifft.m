function [stiffel,nstiffel]=stifft(xy,order,intxi,inteta,intw,curlbast,xmin,ymin,tx,ty,mur,localblendt)

nip=length(intw);
if order==0
esize=3;
else
esize=(order+1)*(order+2);
end
stiffel=zeros(esize);
nstiffel=zeros(esize);
for p=1:nip
%ph=curlbasist(order,intxi(p),inteta(p));
ph=curlbast(:,p);
J=maptri(xy,intxi(p),inteta(p),localblendt);
detJ=det(J);
phxy=ph./detJ;

% PML TENSOR
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
ten=1./(z(2)*z(1));	  
tphxy=ten*phxy;

% non vectorised
% for i=1:esize
% for jj=1:esize
% stiffel(i,jj)=stiffel(i,jj)+(intw(p)*detJ*ten*phxy(i)*phxy(jj)/mur);
% nstiffel(i,jj)=nstiffel(i,jj)+(intw(p)*detJ*phxy(i)*phxy(jj)/mur);
% end
% end

% vectorised
 stiffel=stiffel+intw(p)*detJ*tphxy*phxy'/mur;
 nstiffel=nstiffel+intw(p)*detJ*phxy*phxy'/mur;



end
end
