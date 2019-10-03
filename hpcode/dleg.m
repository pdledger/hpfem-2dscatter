function dln=dleg(x,n)

if n <=0  
dln=0.;
else
a=((-2.*n.*n.*x.*leg(x,n))+(2.*n.*n.*leg(x,n-1)));
b=(2.*n.*(1-(x.^2)));
if a==0 & b==0
dln=0;
elseif a~=0 & b==0
disp('error in gettting dln');
else
dln=a./b;
end
end
