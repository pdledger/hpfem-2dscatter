function [epetq]=splitquad(order2)
   
dx=1./(order2+2);
dy=dx;
epetq=zeros((order2+3)^2,2);
    
np=0;
for i=1:order2+3
   y=0+((i-1)*dy);
   for j=1:order2+3;
      x=0+((j-1)*dx);
      np=np+1;   
      if x >1
        x=1
      elseif x <0
        x=0
      end
      if y >1
        y=1
      elseif y<0
       y=0
      end 
      epetq(np,1)=x;
      epetq(np,2)=y;
   end
end
