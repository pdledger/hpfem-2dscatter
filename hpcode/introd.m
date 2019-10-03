function [sol,unkl,unkh,nunki]=introd(sol,dirval,unkl,unkh,order,ne,nunki)

for i=1:ne
if unkl(i) <0
sol(nunki+abs(unkl(i)))=dirval(abs(unkl(i)));
unkl(i)=nunki+abs(unkl(i));
end
end

for i=1:ne
for p=1:order
if unkh(i,p)<0
sol(nunki+abs(unkh(i,p)))=dirval(abs(unkh(i,p)));
unkh(i,p)=nunki+abs(unkh(i,p));
end
end
end

