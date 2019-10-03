function [npmlb,nzb,nhb,apml,apmlz,apmlh,az,azabs,azh,azpml,ah,ahabs,ahz,ahpml]=extract(ne,order,unkl,unkh,a,aabs,nunk,nunkpmlz,nedgepml)

% to apply the Block Jacobi and Block G-S preconditioners we first need to find
% relevent blocks in the matrix

% use a block containing also the PML functions.
npmlb=nunkpmlz

% construct sparse matrix for zero block
% the global degree of freedom is the same as in this block
apml=a(1:npmlb,1:npmlb);

apmlabs=aabs(1:npmlb,1:npmlb);

% lower order edge block
block=[];
for i=nedgepml+1:ne
  if unkl(i) >0
     block=[block;unkl(i)];
  end
end
% the number of degrees of freedom in this block
nzb=length(block);


% Higher order edge block
block=[];
for i=nedgepml+1:ne
for j=1:order
  if unkh(i,j) >0
     block=[block;unkh(i,j)];
  end
  end
end
% the number of degrees of freedom in this block
nhb=length(block);

if nunk~=npmlb+nzb+nhb
disp('warning incompatable block sizes');
nunk
npmlb,nzb,nhb
end


apmlz=a(1:npmlb,npmlb+1:npmlb+nzb);
apmlh=a(1:npmlb,npmlb+nzb+1:nunk);

az=a(npmlb+1:npmlb+nzb,npmlb+1:npmlb+nzb);
azabs=aabs(1:npmlb:npmlb:nzb,1:npmlb+nzb);
azpml=a(npmlb+1:npmlb+nzb,1:npmlb);
azh=a(npmlb+1:npmlb+nzb,npmlb+nzb+1:nunk);

ah=a(npmlb+nzb+1:nunk,npmlb+nzb+1:nunk);
ahabs=aabs(npmlb+nzb+1:nunk,npmlb+nzb+1:nunk);
ahpml=a(npmlb+nzb+1:nunk,1:npmlb);
ahz=a(npmlb+nzb+1:nunk,npmlb+1:npmlb+nzb);




