function meshsize(Mesht,Meshq,nelemt,nelemq)
coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
intmat=Mesht.Elements;
maxh=0;
for i=1:nelemq
for j=1:4
for k=1:2
xy(j,k)=coord(intmaq(i,j),k);
end
end
l(1)=sqrt((xy(1,1)-xy(2,1))^2+(xy(1,2)-xy(2,2))^2);
l(2)=sqrt((xy(2,1)-xy(3,1))^2+(xy(2,2)-xy(3,2))^2);
l(3)=sqrt((xy(3,1)-xy(4,1))^2+(xy(3,2)-xy(4,2))^2);
l(4)=sqrt((xy(4,1)-xy(1,1))^2+(xy(4,2)-xy(1,2))^2);
avh=(l(1)+l(2)+l(3)+l(4))/4;
if avh > maxh
maxh=avh;
elemno=i;
end
end
if nelemq~=0
disp(['The largest quad element is element',num2str(elemno)]);
disp(['spacing=',num2str(elemno)]);
maxhq=maxh;
end
maxh=0;
for i=1:nelemt
for j=1:3
for k=1:2
xy(j,k)=coord(intmat(i,j),k);
end
end
l(1)=sqrt((xy(1,1)-xy(2,1))^2+(xy(1,2)-xy(2,2))^2);
l(2)=sqrt((xy(2,1)-xy(3,1))^2+(xy(2,2)-xy(3,2))^2);
l(3)=sqrt((xy(3,1)-xy(1,1))^2+(xy(3,2)-xy(1,2))^2);
avh=(l(1)+l(2)+l(3))/3;
if avh > maxh
maxh=avh;
elemno=i;
end
end
disp(['The largest tri element is element',num2str(elemno)]);
disp(['spacing=',num2str(maxh)]);

