function [y]=myproduct(x,a,ahigher,azero,nzb,nhb,maxit,tolh)
t=a*x;
y=blockjac(azero,ahigher,t,nzb,nhb,tolh,maxit);
