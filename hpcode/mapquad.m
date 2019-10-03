function J=mapquad(xy,x,y)

dlam=[-(1-y) -(1-x); (1-y) -x; y x; -y (1-x)];

J=zeros(2);
for i=1:2
for j=1:2
for k=1:4
J(i,j)=J(i,j)+(dlam(k,j)*xy(k,i));
end
end
end
