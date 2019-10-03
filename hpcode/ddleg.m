function ddln=ddleg(x,n)

if n <=1  
ddln=0.;
else
u=((-2.*n.*n.*x.*leg(x,n))+(2.*n.*n.*leg(x,n-1)));
v=(2.*n.*(1-(x.^2)));
du=-2.*n.*n.*leg(x,n)-2.*n.*n.*x.*dleg(x,n)+2.*n.*n.*dleg(x,n-1);
dv=-2.*n.*2.*x;
if v.*du-u.*dv==0 && v.^2==0
ddln=0;
elseif v.*du-u.*dv~=0 && v.^2==0
disp('error in gettting dln');
else
ddln=(v.*du-u.*dv)./v.^2;
end
end
