function [kcc,rc,kccabs]=static1(massel,stiffel,rhsel,order,omega,esize,cont,nstiffel,nmassel)

if esize~=cont
acc=stiffel(1:cont,1:cont)-omega^2*massel(1:cont,1:cont);
aci=stiffel(1:cont,cont+1:esize)-omega^2*massel(1:cont,cont+1:esize);
aic=stiffel(cont+1:esize,1:cont)-omega^2*massel(cont+1:esize,1:cont);
aii=stiffel(cont+1:esize,cont+1:esize)-omega^2*massel(cont+1:esize,cont+1:esize);

lc=rhsel(1:cont);
li=rhsel(cont+1:esize);

kcc=acc-(aci*(aii\aic));
rc=lc-(aci*(aii\li));

% Absolute values for preconditioner
%acc=(stiffel(1:cont,1:cont))+(omega^2*massel(1:cont,1:cont));
%aci=(stiffel(1:cont,cont+1:esize))+(omega^2*massel(1:cont,cont+1:esize));
%aic=(stiffel(cont+1:esize,1:cont))+(omega^2*massel(cont+1:esize,1:cont));
%aii=(stiffel(cont+1:esize,cont+1:esize))+abs(omega^2*massel(cont+1:esize,cont+1:esize));

acc=nstiffel(1:cont,1:cont)+omega^2*nmassel(1:cont,1:cont);
aci=nstiffel(1:cont,cont+1:esize)+omega^2*nmassel(1:cont,cont+1:esize);
aic=nstiffel(cont+1:esize,1:cont)+omega^2*nmassel(cont+1:esize,1:cont);
aii=nstiffel(cont+1:esize,cont+1:esize)+omega^2*nmassel(cont+1:esize,cont+1:esize);

kccabs=acc-(aci*(aii\aic));

else
kcc=stiffel-omega^2.*massel;
kccabs=nstiffel+omega^2.*nmassel;
rc=rhsel;
end
