function [basisx,basisy,curlbasis,rbasisx,rbasisy]=myeval(order,x,w)

nip=length(x);

basisx=[];
basisy=[];
curlbasis=[];
for p=1:nip
for q=1:nip
ph=basisq(order,x(p),x(q));
basisx=[basisx, ph(:,1)];
basisy=[basisy, ph(:,2)];
cph=curlbasisq(order,x(p),x(q));
curlbasis=[curlbasis, cph];
end
end

% 4 bc edge possibilities
% Edge 1 y=0
rbasisx=[];
rbasisy=[];
for p=1:nip
ph=basisq(order,x(p),0);
rbasisx=[rbasisx, ph(:,1)];
rbasisy=[rbasisy, ph(:,2)];
end
% Edge 2 y=1
for p=1:nip
ph=basisq(order,x(p),1);
rbasisx=[rbasisx, ph(:,1)];
rbasisy=[rbasisy, ph(:,2)];
end
% Edge 3 x=1
for p=1:nip
ph=basisq(order,1,x(p));
rbasisx=[rbasisx, ph(:,1)];
rbasisy=[rbasisy, ph(:,2)];
end
% Edge 4 x=0
for p=1:nip
ph=basisq(order,0,x(p));
rbasisx=[rbasisx, ph(:,1)];
rbasisy=[rbasisy, ph(:,2)];
end


