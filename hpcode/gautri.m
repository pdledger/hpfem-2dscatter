function [intxi,inteta,intw]=gautri(n)

% first compute the points for the integral int(-1,1) (1+x) f(x) dx
% use gaujac with alf=0,bet=1
alf=0;
bet=1;
[xj,wj]=gaujac(alf,bet,n);
% next transform the intergal to int(0,1) (1-x) dx
u=(1-xj)./2;
b=wj./4;

% use gaujac with alf=0,bet=0
alf=0;
bet=0;
[xj,wj]=gaujac(alf,bet,n);
% next transform the intergal to int(0,1)  dx
v=(1+xj)./2;
c=wj./2;

% compute weights and locations
p=0;
for i=1:n
for j=1:n
  p=p+1;
  intxi(p)=u(i);
  inteta(p)=v(j)*(1-u(i));
  intw(p)=b(i)*c(j);
end
end
nip=p;

% Area of our reference triangle is sqrt(3). readjust weights so that they integrate a
%constant (previous triangle has area 1/2)
intw=intw.*2*sqrt(3);
% compute locations on master element
xy=[1. 0.;
    0. sqrt(3.);
    -1. 0.];
for i=1:nip
   xi=0.;
   eta=0.;
   l(1)=1-intxi(i)-inteta(i);
   l(2)=intxi(i);
   l(3)=inteta(i);
   for j=1:3
      xi=xi+(xy(j,1)*l(j));
      eta=eta+(xy(j,2)*l(j));
   end
   intxi(i)=xi;
   inteta(i)=eta;
end
