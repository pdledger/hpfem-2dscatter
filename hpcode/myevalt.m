function [basisxt,basisyt,curlbast,rbasisxt,rbasisyt]=myevalt(...
          order,intxi,inteta,xt)

nipt=length(intxi);
nip=length(xt);

basisxt=[];
basisyt=[];
curlbast=[];
for p=1:nipt
ph=basist(order,intxi(p),inteta(p));
basisxt=[basisxt, ph(:,1)];
basisyt=[basisyt, ph(:,2)];
cph=curlbasist(order,intxi(p),inteta(p));
curlbast=[curlbast, cph];
end

% 3 bc edge possibilities
% Edge 1 y=0
rbasisxt=[];
rbasisyt=[];
for p=1:nip
% edge 1
eta=sqrt(3.)*0.5*(xt(p)+1);
xi=(-1+xt(p))*(-0.5);
ph=basist(order,xi,eta);
rbasisxt=[rbasisxt, ph(:,1)];
rbasisyt=[rbasisyt, ph(:,2)];
end
% Edge 2 y=1
for p=1:nip
eta=sqrt(3.)*(-0.5)*(xt(p)-1);
xi=-(1+xt(p))*0.5;
ph=basist(order,xi,eta);
rbasisxt=[rbasisxt, ph(:,1)];
rbasisyt=[rbasisyt, ph(:,2)];
end
% Edge 3 x=1
for p=1:nip
eta=0.;
xi=xt(p);
ph=basist(order,xi,eta);
rbasisxt=[rbasisxt, ph(:,1)];
rbasisyt=[rbasisyt, ph(:,2)];
end



