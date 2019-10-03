syms x y l1 l2 l3 l4 s1 s2 s3 s4 ph1 ph2 ph3 ph4;
l1=(1-x)*(1-y);
l2=x*(1-y);
l3=x*y;
l4=(1-x)*y;

s1=(1-x)+(1-y);
s2=x+1-y;
s3=x+y;
s4=(1-x)+y;

ph1=simple([diff(s2-s1,x); diff(s2-s1,y)]*(l1+l2))

ph2=simple([diff(s4-s3,x); diff(s4-s3,y)]*(l4+l3))

ph3=simple([diff(s3-s2,x); diff(s3-s2,y)]*(l2+l3))

ph4=simple([diff(s1-s4,x); diff(s1-s4,y)]*(l4+l1))

eps1=simple(s2-s1)

sig1=simple(l1+l2)

eps2=simple(s4-s3)

sig2=simple(l4+l3)

eps3=simple(s3-s2)

sig3=simple(l3+l2)

eps4=simple(s1-s4)

sig4=simple(l1+l4)
