syms x y l1 l2 l3 l4 s1 s2 s3 s4 ph1 ph2 ph3 ph4;
l1=(sqrt(3)+sqrt(3)*x-y)/(2*sqrt(3));
l2=y/sqrt(3);
l3=(sqrt(3)-sqrt(3)*x-y)/(2*sqrt(3));

s1=simple(l2-l1)
s2=simple(l3-l2)
s3=simple(l1-l3)

t1=simple(l1+l2)
t2=simple(l3+l2)
t3=simple(l1+l3)

ph1=simple([diff(l2,x); diff(l2,y)]*(l1)-[diff(l1,x); diff(l1,y)]*(l2))

ph2=simple([diff(l3,x); diff(l3,y)]*(l2)-[diff(l2,x); diff(l2,y)]*(l3))

ph3=simple([diff(l1,x); diff(l1,y)]*(l3)-[diff(l3,x); diff(l3,y)]*(l1))

