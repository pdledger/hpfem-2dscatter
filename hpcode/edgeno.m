function [ne,edges,glob,dir]=edgeno(Mesht,Meshq,bsido,nelemt,nelemq,nboun,probdata)

% This code could be improved by using linked lists. Basic version supplied here
% only!

intmaq=Meshq.Elements;
intmat=Mesht.Elements;
coord=Mesht.Coordinates;
fnum=Mesht.Fnum;
[npoin dum]=size(coord);
glob=zeros(nelemq+nelemt,4);
dir=zeros(nelemq+nelemt,4);
help1=zeros(npoin,1);

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
      nmax=max([n1 n2]);

      for kk=1:help1(nmax)
      k=help2(nmax,kk);      
%         for k=1:ne
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

% define help arrays to reduce cost of search
if n1 > n2
help1(n1)=help1(n1)+1;
help2(n1,help1(n1))=ne;
else
help1(n2)=help1(n2)+1;
help2(n2,help1(n2))=ne;
end
	 
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

% correct for multiple materials
% helpedge=zeros(ne,1);
% helpedge2=zeros(ne,2);
% for i=1:nelemt
% for j=1:3
% helpedge(glob(i,j))=helpedge(glob(i,j))+1;
% helpedge2(glob(i,j),helpedge(glob(i,j)))=fnum(i);
% end
% end
% for i=1:ne
% if helpedge(i)==2 & helpedge2(i,1)~=helpedge2(i,2)
% % This is an edge on the interface
% edges(i,3)=1;
% end
% end
% correct for multiple materials
helpedge=zeros(ne,1);
helpedge2=zeros(ne,2);
for i=1:nelemt
for j=1:3
helpedge(glob(i,j))=helpedge(glob(i,j))+1;
helpedge2(glob(i,j),helpedge(glob(i,j)))=fnum(i);
end
end
for i=1:ne
if helpedge(i)==2 & helpedge2(i,1)~=helpedge2(i,2)
% these vairables only defined for multiple materials.
nint=probdata.nint;
bcint=probdata.bcint;
    found=-100000;
    for j=1:nint(helpedge2(i,1))
       for k=1: nint(helpedge2(i,2)) 
           if bcint(helpedge2(i,1),j)==bcint(helpedge2(i,2),k)
                found=bcint(helpedge2(i,1),j);
           end
       end
    end
    if found ~=-100000;
        % This is a BC interface
        edges(i,3)=found;
    else
        disp('error')
    end
end
end
