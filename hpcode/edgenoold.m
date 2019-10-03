function [ne,edges,glob,dir]=edgeno(Mesht,Meshq,bsido,nelemt,nelemq,nboun)

% This code could be improved by using linked lists. Basic version supplied here
% only!

intmaq=Meshq.Elements;
intmat=Mesht.Elements;
glob=zeros(nelemq+nelemt,4);
dir=zeros(nelemq+nelemt,4);

ne=0;
for i=1:nelemq+nelemt 

   if i > nelemq
      temp=[intmat(i-nelemq,1) intmat(i-nelemq,2);
            intmat(i-nelemq,2) intmat(i-nelemq,3);
            intmat(i-nelemq,3) intmat(i-nelemq,1)];
      nepel=3;
   else
      temp=[intmaq(i,1) intmaq(i,2);
            intmaq(i,3) intmaq(i,4);
            intmaq(i,2) intmaq(i,3);
            intmaq(i,4) intmaq(i,1)];  
      nepel=4;
   end

   for j=1:nepel      
      n1=temp(j,1);
      n2=temp(j,2);
      flag=0;
      bctyp=0;
      dir1=0;

% check to see if edge is already defined
      if ne~=0
         for k=1:ne
         i1=edges(k,1);
         i2=edges(k,2);
         if (n1==i1 & n2==i2) | (n1==i2 & n2== i1)
           flag=k;
           if n1==i2 
           dir1=1;
           end
         end
         end
      end

% check for a BC value
      if flag==0
         for k=1:nboun
         i1=bsido(k,1);
         i2=bsido(k,2);
         if (n1==i1 & n2==i2) | (n1==i2 & n2==i1)
            bctyp=bsido(k,4);
         end
	 end
      end

      if flag==0
         ne=ne+1;
         edges(ne,1)=n1;
         edges(ne,2)=n2;
         edges(ne,3)=bctyp;
         glob(i,j)=ne;
         dir(i,j)=1.;
         if dir1==1
	    dir(i,j)=-1.;
	 end
      else
      dir(i,j)=1.;
      if dir1==1
        dir(i,j)=-1.;
      end
      glob(i,j)=flag;
      end

end
end      
disp(['We have found',num2str(ne),'edges'])

