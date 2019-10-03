function [bsido,p] = my_data(p,t,cnect,bcflags,node)
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
[e,i,j] = unique(sort(e,2),'rows');                         % Unique edges
vec     = vec(j);                                           % Unique edge numbers
eINt    = [vec(vect), vec(vect+numt), vec(vect+2*numt)];    % Unique edges in each triangle



% DETERMINE EDGE TO TRIANGLE CONNECTIVITY

% Each row has two entries corresponding to the triangle numbers
% associated with each edge. Boundary edges have one entry = 0.
nume = size(e,1);
e2t  = repmat(0,nume,2);
ndx  = repmat(1,nume,1);
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
v = repmat(0,max(j),1);
for k = 1:length(j)
    jk = j(k); v(jk) = v(jk)+1;
end
maxN = max(v);

n2e = repmat(0,numn,maxN+1);
ndx = repmat(1,numn,1);
for k = 1:nume
    % End nodes
    n1 = e(k,1); n2 = e(k,2);
    % Connectivity
    n2e(n1,ndx(n1)) = k; ndx(n1) = ndx(n1)+1;
    n2e(n2,ndx(n2)) = k; ndx(n2) = ndx(n2)+1;
end


% DETERMINE NODE TO NODE CONNECTIVITY

n2n = repmat(0,numn,maxN+1);
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
bnd     = false(numn,1);
bnd(bn) = true;                 % True for boundary nodes

[nbe dum]=size(bn);

disp(['we have',num2str(nbe),'boundary edges'])

% create coordinates of all boundary nodes
xbn=p(bno,1);
ybn=p(bno,2);

% check boundary flag of each boundary point
[nbs dum]=size(cnect);
bcn=zeros(numn,2);
TOL=1e-6;
help=zeros(numn,1);
for i=1:length(xbn)

    for j=1:nbs
        xmin=min(node(cnect(j,1),1),node(cnect(j,2),1));
        xmax=max(node(cnect(j,1),1),node(cnect(j,2),1));

        ymin=min(node(cnect(j,1),2),node(cnect(j,2),2));
        ymax=max(node(cnect(j,1),2),node(cnect(j,2),2));
        
        x1=node(cnect(j,1),1);
        y1=node(cnect(j,1),2);
        x2=node(cnect(j,2),1);
        y2=node(cnect(j,2),2);

        if (xmax-xmin)^2+(ymin-ymax)^2 == 0
            disp('boundary segment of zero length!');
        end

        out=0;
        if xbn(i) >= xmin-TOL && xbn(i) <= xmax+TOL && ybn(i) >= ymin-TOL && ybn(i)<=ymax+TOL
            % this point lies within the min and max points
            %disp('found a point')

            if abs(ymin-ymax) < TOL
                if abs(ybn(i)-ymax) < TOL
                    out=1;
                end
            elseif abs(xmin-xmax) < TOL
                if abs(xbn(i)-xmax) < TOL
                    out=1;
                end
            else
                
                m=(y2-y1)/(x2-x1);
                c=y2-m*x2;
 %               m=(ymax-ymin)/(xmax-xmin);
 %               c=ymax-m*xmax;
                if abs(ybn(i)-(m*xbn(i)+c)) < TOL
                    out=1;
                end
            end
        end
        if out==1;
            help(bno(i))=help(bno(i))+1;
            if help(bno(i)) > 2
                disp('more than 2 edges meet at a node!!')
            end
            bcn(bno(i),help(bno(i)))=bcflags(j);
        end
    end

    if help(bno(i))==0
        xbn(i),ybn(i)
        node
        pause
    end

end



% plot out boundary conditions
figure
hold on
for i=1:length(xbn)
    if bcn(bno(i),1) == 0
        % no boundary condition allocated
        plot(xbn(i),ybn(i),'bx')
    elseif bcn(bno(i),1) == 1
        plot(xbn(i),ybn(i),'ro')
    elseif bcn(bno(i),1) == 2
        plot(xbn(i),ybn(i),'yo')
    elseif bcn(bno(i),1) == 3
        plot(xbn(i),ybn(i),'go')
    elseif bcn(bno(i),1) == 4
        plot(xbn(i),ybn(i),'ko')
    end
end
hold off


for i=1:nbe
    bsido(i,1)=bn(i,1);
    bsido(i,2)=bn(i,2);

    if help(bn(i,1))==1;
        bcleft1=bcn(bn(i,1),1);
        bcleft2=0;
    else
        bcleft1=bcn(bn(i,1),1);
        bcleft2=bcn(bn(i,1),2);
    end

    if help(bn(i,2))==1;
        bcright1=bcn(bn(i,2),1);
        bcright2=0;
    else
        bcright1=bcn(bn(i,2),1);
        bcright2=bcn(bn(i,2),2);
    end

    if bcleft1==bcright1
        bsido(i,3)=0;
        bsido(i,4)=bcleft1;
        bsido(i,5)=bcleft1;
    elseif bcleft2~=0 & bcright2==0
        bsido(i,3)=0;
        bsido(i,4)=bcright1;
        bsido(i,5)=bcright1;
    elseif bcleft2==0 & bcright2~=0
        bsido(i,3)=0;
        bsido(i,4)=bcleft1;
        bsido(i,5)=bcleft1;
    else
        disp('not possible')
    end

end


