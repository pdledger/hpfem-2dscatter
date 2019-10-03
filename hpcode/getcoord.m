function [x,y]=getcoord(xy,xi,et)
la=[(1-xi)*(1-et); xi*(1-et);xi*et; (1-xi)*et];
x=0;
y=0;
for j=1:4
  x=x+(xy(j,1)*la(j));
  y=y+(xy(j,2)*la(j));
end

