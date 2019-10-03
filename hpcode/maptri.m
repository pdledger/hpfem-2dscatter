function J=maptri(xy,epsi,eta,localblend)

%dlam=(1/sqrt(3))*[sqrt(3)/2 -1/2 ; 0 1; -sqrt(3)/2 -1/2];

%J=zeros(2);
%for i=1:2
%for j=1:2
%for k=1:3
%J(i,j)=J(i,j)+(dlam(k,j)*xy(k,i));
%end
%end
%end


l(1)=(1./(2*sqrt(3.)))*(sqrt(3.)+(sqrt(3.)*epsi)-eta);
l(2)=eta/sqrt(3.);
l(3)=(1./(2*sqrt(3.)))*(sqrt(3.)-(sqrt(3.)*epsi)-eta);
      
%Derivitives of area coordinates      
dl(1,1)=0.5;
dl(1,2)=-1.*(1./(2*sqrt(3.)));
dl(2,1)=0.;
dl(2,2)=1./sqrt(3.);
dl(3,1)=-0.5;
dl(3,2)=-(1./(2*sqrt(3.)));
      
dxdl(1)=xy(1,1)+(l(2)*... 
gx(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2)))+...          
(l(1)*l(2)*dgx(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2))*(-1.))+...
(l(3)*gx(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2)))+...
(l(3)*l(1)*dgx(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2))*(1.));

dxdl(2)=xy(2,1)+(l(1)*gx(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2)))+...
(l(1)*l(2)*dgx(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2))*(1.))+...
(l(3)*gx(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2)))+...
(l(3)*l(2)*dgx(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2))*(-1.));

dxdl(3)=xy(3,1)+(l(2)* ...
gx(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2)))+...
(l(2)*l(3)*dgx(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2))*(1.))+...
(l(1)*gx(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2)))+...
(l(3)*l(1)*dgx(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2))*(-1.));
 

dydl(1)=xy(1,2)+(l(2)* ...
gy(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2)))+ ...
(l(1)*l(2)*dgy(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2))*(-1.))+...
(l(3)*gy(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2)))+ ...
(l(3)*l(1)*dgy(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2))*(1.));
     
     
dydl(2)=xy(2,2)+(l(1)* ...
gy(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2)))+ ...
(l(1)*l(2)*dgy(l(2)-l(1),localblend(1,1),localblend(1,2),xy(1,1),xy(2,1),xy(1,2),xy(2,2))*(1.))+...
(l(3)*gy(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2)))+      ...
(l(3)*l(2)*dgy(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2))*(-1.));
     
dydl(3)=xy(3,2)+(l(2)*...
gy(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2)))+  ...
(l(2)*l(3)*dgy(l(3)-l(2),localblend(2,1),localblend(2,2),xy(2,1),xy(3,1),xy(2,2),xy(3,2))*(1.))+  ...
(l(1)*gy(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2)))+ ...
(l(3)*l(1)*dgy(l(1)-l(3),localblend(3,1),localblend(3,2),xy(3,1),xy(1,1),xy(3,2),xy(1,2))*(-1.));

 
%     jac(1,1)=dx/depsi
J(1,1)=0.;
      for j=1:3
J(1,1)=J(1,1)+(dxdl(j)*dl(j,1));
      end
      
%     jac(1,2)=dx/deta      
J(1,2)=0.;
      for j=1:3
J(1,2)=J(1,2)+(dxdl(j)*dl(j,2));
      end
      
      
%     jac(2,1)=dy/depsi      
J(2,1)=0.;
      for j=1:3
J(2,1)=J(2,1)+(dydl(j)*dl(j,1));
      end
      
%     jac(2,2)=dy/deta      
J(2,2)=0.;
      for j=1:3
J(2,2)=J(2,2)+(dydl(j)*dl(j,2));
      end

if det(J) < 0 
dlam=(1/sqrt(3))*[sqrt(3)/2 -1/2 ; 0 1; -sqrt(3)/2 -1/2];

J=zeros(2);
for i=1:2
for j=1:2
for k=1:3
J(i,j)=J(i,j)+(dlam(k,j)*xy(k,i));
end
end
end
disp('linear geometry used');
end
      
      
