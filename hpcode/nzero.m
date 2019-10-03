function [nz,map,irn,icn]=nzero(unkl,unkh,unkqi,unkti,nelemq,nelemt,order,glob,nunki)

%create the connectivity table for quadrilaterals
connectq=[];
for i=1:nelemq
data=[];
% Low order block
for j=1:4
data=[data unkl(glob(i,j))];
end
% High order block
for j=1:4
for p=1:order
data=[data unkh(glob(i,j),p)];
end
end
connectq=[connectq; data ];
end
ldataq=4*(order+1);
disp('Found all quadrilateral connections')


%create the connectivity table for triangles
connectt=[];
for i=1:nelemt
data=[];
%% Low order block
for j=1:3
data=[data unkl(glob(i+nelemq,j))];
end
% High order block
for j=1:3
for p=1:order
data=[data unkh(glob(i+nelemq,j),p)];
end
end
% Interior block
%for j=1:order*(order-1)+order-1
%data=[data unkti(i,j)];
%end
connectt=[connectt; data ];
end
ldatat=3*(order+1);
%if order==0
%ldatat=3;
%else
%ldatat=(order+1)*(order+2);
%end
disp('Found all triangular connections')

% to here !!!


% create the structure ready for building the nonzero enteries
help1=zeros(nunki,1); % help(i)=Number of Unknowns connected to Unknown i
help2=zeros(nunki,ldataq); % help(i,j)=List of Unknowns connected to Unknown i
help3=zeros(nunki,1);
for i=1:nelemq
   for j=1:ldataq
      Row=connectq(i,j);
      if Row > 0 
      
%      transfer enteries in to help3
       for q=1:help1(Row)
       help3(help2(Row,q))=1;
       end
      
         for k=1:ldataq
            Col=connectq(i,k);
            if Col > 0
% check to see if Columne Node is already stored
%               flag=0;
%               if( help1(Row) > 0)
%                  for p=1:help1(Row)
%                     if (help2(Row,p) == Col) 
%                        flag=1;
%                     end
%                  end
%               end
%               if(flag==0)
		if(help3(Col)~=1)
                  help1(Row)=help1(Row)+1;
                  help2(Row,help1(Row))=Col;
               end
	    end
         end
	 %        reintialise help array
         for q=1:help1(Row)
         help3(help2(Row,q))=0;
         end
	

      end
   end
end

for i=1:nelemt
   for j=1:ldatat
      Row=connectt(i,j);
      if Row > 0 
         for k=1:ldatat
            Col=connectt(i,k);
            if Col > 0
% check to see if Columne Node is already stored
               flag=0;
               if( help1(Row) ~= 0)
                  for p=1:help1(Row)
                     if (help2(Row,p) == Col) 
                        flag=1;
                     end
                  end
               end
               if(flag==0)
                  help1(Row)=help1(Row)+1;
                  help2(Row,help1(Row))=Col;
               end
	    end
         end
      end
   end
end

disp('finished part1')
% create the list of Nonzero enteries in the matrix
nz=0;
map=zeros(nunki+1,1); %map(i)=First non-zero entry on row i
for i=1:nunki
%store the nonzero entry on which each row begins
   map(i)=nz+1;
   for j=1:help1(i)
      nz=nz+1;
      irn(nz)=i;
      icn(nz)=help2(i,j);
   end
end
map(nunki+1)=nz+1;

disp(['We have found'])
nz
clear help1
clear help2
clear connectq
clear connectt
