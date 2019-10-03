function [bsido,p] = my_data_axi(p,t,cnect,bcflags,node)
%scatterer,bcscatter)

% Build data structures used in the solver.
%
% This function builds the mesh parameters and connectivity that will be
% used in tvd_rk2.m. Building this initially saves HEAPS of CPU time later
% in tvd_rk2.m, but of course, my meshes are fixed for the integration.
%
% Darren Engwirda - 2005-2006.
%
% Naver2d is Copyright (C) 2005-2006 Darren Engwirda. See "copyright.m" for
% full details.

numn = size(p,1);
numt = size(t,1);
vect = 1:numt;

% DETERMINE UNIQUE EDGES IN MESH

e       = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];             % Edges - not unique
vec     = (1:size(e,1))';                                   % List of edge numbers
[e,~,j] = unique(sort(e,2),'rows');                         % Unique edges
vec     = vec(j);                                           % Unique edge numbers
eINt    = [vec(vect), vec(vect+numt), vec(vect+2*numt)];    % Unique edges in each triangle


% DETERMINE EDGE TO TRIANGLE CONNECTIVITY

% Each row has two entries corresponding to the triangle numbers
% associated with each edge. Boundary edges have one entry = 0.
nume = size(e,1);
% e2t  = repmat(0,nume,2);
% ndx  = repmat(1,nume,1);
e2t  = zeros(nume,2);
ndx  = ones(nume,1);
for k = 1:numt
    % Edge in kth triangle
    e1 = eINt(k,1); e2 = eINt(k,2); e3 = eINt(k,3);
    % Edge 1
    e2t(e1,ndx(e1)) = k; ndx(e1) = ndx(e1)+1;
    % Edge 2
    e2t(e2,ndx(e2)) = k; ndx(e2) = ndx(e2)+1;
    % Edge 3
    e2t(e3,ndx(e3)) = k; ndx(e3) = ndx(e3)+1;
end



% DETERMINE NODE TO EDGE CONNECTIVITY

% Determine maximum neighbours
j = e(:);
v = zeros(max(j),1);
for k = 1:length(j)
    jk = j(k); v(jk) = v(jk)+1;
end
maxN = max(v);

n2e = zeros(numn,maxN+1);
ndx = ones(numn,1);
for k = 1:nume
    % End nodes
    n1 = e(k,1); n2 = e(k,2);
    % Connectivity
    n2e(n1,ndx(n1)) = k; ndx(n1) = ndx(n1)+1;
    n2e(n2,ndx(n2)) = k; ndx(n2) = ndx(n2)+1;
end


% DETERMINE NODE TO NODE CONNECTIVITY

n2n = zeros(numn,maxN+1);
for k = 1:numn
    next = 1; m = 1;
    while n2e(k,m)>0
        if e(n2e(k,m),1)==k
            n2n(k,next) = e(n2e(k,m),2); next = next+1;
        else
            n2n(k,next) = e(n2e(k,m),1); next = next+1;
        end
        m = m+1;
    end
end


% FLAG BOUNDARY ELEMENTS
vec     = (1:nume)';            % Edge list
be      = vec(~all(e2t,2));     % Boundary edges
bn      = e(be,:);
bno     = unique(bn(:));        % Boundary nodes
% bnd     = false(numn,1);
% bnd(bn) = true;                 % True for boundary nodes

nbe=size(bn,1);
disp(['We have ',num2str(nbe),' boundary edges'])

% create coordinates of all boundary nodes
xbn=p(bno,1);
ybn=p(bno,2);

% check boundary flag of each boundary point
nbs=size(cnect,1); 
bce=zeros(size(bn)); 
TOL=1e-6; 
TOL2=1e-10;
TOLn= 1e-12;

for i=1:nbe
    % Find vector equation of the line on the mesh (v = a*i+b*j+0*k)
    vecm = [p(bn(i,2),1) - p(bn(i,1),1); p(bn(i,2),2) - p(bn(i,1),2)];
    vecm = vecm/norm(vecm);

        x1=p(bn(i,1),1);             y1=p(bn(i,1),2);
	x2=p(bn(i,2),1);             y2=p(bn(i,2),2);
        xmaxg=max(x1,x2);
	xming=min(x1,x2);
	ymaxg=max(y1,y2);
	yming=min(y1,y2);
	

    for j=1:nbs
        x1=node(cnect(j,1),1);        y1=node(cnect(j,1),2);
        x2=node(cnect(j,2),1);        y2=node(cnect(j,2),2);
        xmaxbs=max(x1,x2);
	xminbs=min(x1,x2);
	ymaxbs=max(y1,y2);
	yminbs=min(y1,y2);
	
        % Find vector equation of the line on the geometry
        vecg = [x2-x1; y2-y1];
        vecg=vecg/norm(vecg);
        
        angle = acos(vecg'*vecm);
        
	
	if xming >= xminbs-TOL2 & xmaxg <= xmaxbs+TOL2 & yming >= yminbs-TOL2 & ymaxg <= ymaxbs+TOL2
                    
        	if abs(angle) < TOL || abs(angle) < pi-TOL || abs(angle) < pi+TOL 
                
            	% Now find the distance between the two parallel lines
            	% Find the difference vector
%            	vecd = [x1-xbn(i); y1-ybn(i)]; 
            	vecd = [x1-p(bn(i,1),1); y1-p(bn(i,1),2)]; 
             	if norm(vecd) < TOLn
            	  % X = [xbn(i+1), ybn(i+1); x2, y2];
%            	  X = [xbn(i), ybn(i); x1, y1];
            	  X = [p(bn(i,1),1), p(bn(i,1),2); x1, y1];

%             	  d = octave_pdist(X,'euclidean'); 
                d = pdist(X,'euclidean'); 
            	else
               	  vecd=vecd/norm(vecd);
            	  %  vpro =(vecm'*vecd)*vecm;
            	  vpro =(vecg'*vecd)*vecg;
            	  d = norm(vecd-vpro);
            	end
        
       	        if d < TOL
                  % Lines coincide
                  bce(i,:)=bcflags(j); 
		  %xming, xminbs,xmaxg,xmaxbs,
		  %yming, yminbs,ymaxg,ymaxbs
		  %bcflags(j)
		  %pause
		  
                end
            end
    end
    end    
end

bsido = zeros(nbe,5);
bsido(:,1:2)=bn;
bsido(:,4:5)=bce;

%plot out boundary conditions
figure
hold on
for i=1:nbe
    for j=1:2
    if bsido(i,4) == 0
        % no boundary condition allocated
        plot(p(bn(i,j),1),p(bn(i,j),2),'bx')
    elseif bsido(i,4) == 1
        plot(p(bn(i,j),1),p(bn(i,j),2),'ro')
    elseif bsido(i,4) == 2
        plot(p(bn(i,j),1),p(bn(i,j),2),'yo')
    elseif bsido(i,4) == 3
        plot(p(bn(i,j),1),p(bn(i,j),2),'go')
    elseif bsido(i,4) == 4
        plot(p(bn(i,j),1),p(bn(i,j),2),'ko')
    end
    end
end


% for i=1:length(xbn)
%    if bsido(i,4) == 0
%        % no boundary condition allocated
%        plot(xbn(i),ybn(i),'bx')
%    elseif bsido(i,4) == 1
%        plot(xbn(i),ybn(i),'ro')
%    elseif bsido(i,4) == 2
%        plot(xbn(i),ybn(i),'yo')
%    elseif bsido(i,4) == 3
%        plot(xbn(i),ybn(i),'go')
%    elseif bsido(i,4) == 4
%        plot(xbn(i),ybn(i),'ko')
%    end
% end
hold off
% legend('No BC','DirX-DirY','DirX - NeuY', 'NeuX - DirY', 'NeuX-NeuY')
