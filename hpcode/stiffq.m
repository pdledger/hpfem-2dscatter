function [stiffel,nstiffel]=stiffq(xy,order,x,w,curlbasis,xmin,ymin,tx,ty,localblend)

nip=length(x);
esize=4*(order+1)+2*order*(order+1);
stiffel=zeros(esize);
nstiffel=zeros(esize);
for p=1:nip
for q=1:nip
ph=curlbasis(:,q+(nip*(p-1)));
J=mapquad(xy,x(p),x(q),localblend);
detJ=det(J);
phxy=ph./detJ;

% PML TENSOR
[xp,yp]=getcoord(xy,x(p),x(q),localblend);
c=-log(eps);
if abs(xp) >= xmin
%z(1)=1-j*(10*((abs(xp)-xmin)/0.5).^4*(xp*sign(xp))+...
%          ((abs(xp)-xmin)/0.5)^5);
z(1)=1-5*j*c*(abs(xp)-xmin)^4/tx^5*xp*sign(xp)-j*c*(abs(xp)-xmin)^5/tx^5;
else
z(1)=1;
end
if abs(yp) >= ymin
%z(2)=1-j*(10*((abs(yp)-ymin)/0.5).^4*(yp*sign(yp))+...
%          ((abs(yp)-ymin)/0.5)^5);
z(2)=1-5*j*c*(abs(yp)-ymin)^4/ty^5*yp*sign(yp)-j*c*(abs(xp)-ymin)^5/ty^5;
else
z(2)=1;
end
ten(1)=1./(z(2).*z(1));	  


for i=1:esize
for jj=1:esize
stiffel(i,jj)=stiffel(i,jj)+(w(p)*w(q)*detJ*ten(1)*phxy(i)*phxy(jj));
nstiffel(i,jj)=nstiffel(i,jj)+(w(p)*w(q)*detJ*phxy(i)*phxy(jj));
end
end

end
end
