function mycontourz(nelemt,nelemq,orders,Meshst,Meshsq,order,splitq,unkl,unkh...
    ,unkqi,unkti,dir,glob,Meshq,Mesht,sol,splitt,coefft,coeffq,omega)

coordn=Meshsq.Coordinates;
[npoinn,i]=size(coordn);
intmaq=Meshq.Elements;
intmat=Mesht.Elements;
coord=Meshq.Coordinates;

% Quantities for triangular splittings
nintpt=0;
for i=1:orders
    nintpt=nintpt+i;
end
%  number of elements per triangle
nelpt=4;
for i=1:orders
    nelpt=nelpt+3+(2*i);
end
% number of points per triangle
nptri=6;
for i=1:orders
    nptri=nptri+(3+i);
end
epetq=splitquad(orders);
[epett,layer]=splittri(orders,nptri);

phxstore=[];
phystore=[];
cphstore=[];
for q=1:(orders+3)^2
    xi=epetq(q,1);
    eta=epetq(q,2);
    ph=basisq(order,xi,eta);
    phxstore=[phxstore, ph(:,1)];
    phystore=[phystore, ph(:,2)];
    cph=curlbasisq(order,xi,eta);
    cphstore=[cphstore, cph(:)];
end
phxtstore=[];
phytstore=[];
cphtstore=[];
for q=1:nptri
    xi=0.99*epett(q,1);
    eta=0.99*epett(q,2);
    ph=basist(order,xi,eta);
    phxtstore=[phxtstore, ph(:,1)];
    phytstore=[phytstore, ph(:,2)];
    cph=curlbasist(order,xi,eta);
    cphtstore=[cphtstore, cph(:)];
end

out=zeros(npoinn,2);
help=zeros(npoinn,1);

%     loop over quad elements...
for j=1:nelemq
    esize=4*(order+1)+2*order*(order+1);
    for q=1:(orders+3)^2
        for p=1:4
            for pp=1:2
                xy(p,pp)=coord(intmaq(j,p),pp);
            end
        end
        
        for p=1:4
            for pp=1:2
                localblend(p,pp)=coeff(j,p,pp);
            end
        end
        
        xi=epetq(q,1);
        eta=epetq(q,2);
        J=mapquad(xy,xi,eta,localblend);
        detJ=det(J);
        Jinvt=(inv(J))';
        phx=phxstore(:,q);
        phy=phystore(:,q);
        ph=[phx phy];
        phxy=Jinvt*ph';
        phxy=phxy';
        cph=cphstore(:,q);
        cphxy=cph./detJ;
        % lowest order block
        for l=1:4
            lunk(l)=unkl(glob(j,l));
            ldr(l)=dir(j,l);
        end
        % Higher order block
        if order> 0
            for l=1:4
                for p=1:order
                    lunk(l+4*p)=unkh(glob(j,l),p);
                    ldr(l+4*p)=dir(j,l)^(p-1);
                end
            end
            %Interior block
            for l=1:2*order*(order+1)
                lunk(4*(order+1)+l)=unkqi(j,l);
                ldr(4*(order+1)+l)=1;
            end
        end
        etil=zeros(2,1);
        rcetil=0;
        icetil=0;
        
        for k=1:esize
            if lunk(k) ~=0
                etil(1)=etil(1)+(phxy(k,1)'.*(sol(lunk(k))).*ldr(k));
                etil(2)=etil(2)+(phxy(k,2)'.*(sol(lunk(k))).*ldr(k));
                
            end
        end
        
        out(splitq(j,q),1)=out(splitq(j,q),1)+etil(1);
        out(splitq(j,q),2)=out(splitq(j,q),2)+etil(2);
        help(splitq(j,q))=help(splitq(j,q))+1;
    end
end

%     loop over triangular elements...
for j=1:nelemt
    esize=3;
    if order>0
        esize=(order+1)*(order+2);
    end
    for q=1:nptri
        for p=1:3
            for pp=1:2
                xyt(p,pp)=coord(intmat(j,p),pp);
            end
        end
        
        for p=1:3
            for pp=1:2
                localblendt(p,pp)=coefft(j,p,pp);
            end
        end
        
        xi=epett(q,1);
        eta=epett(q,2);
        J=maptri(xyt,xi,eta,localblendt);
        detJ=det(J);
        Jinvt=(inv(J))';
        phx=phxtstore(:,q);
        phy=phytstore(:,q);
        ph=[phx phy];
        phxy=Jinvt*ph';
        phxy=phxy';
        cph=cphtstore(:,q);
        cphxy=cph./detJ;
        
        % lowest order block
        for l=1:3
            lunk(l)=unkl(glob(j+nelemq,l));
            ldr(l)=dir(j+nelemq,l);
        end
        % Higher order block
        if order> 0
            for l=1:3
                for p=1:order
                    lunk(l+3*p)=unkh(glob(j+nelemq,l),p);
                    ldr(l+3*p)=dir(j+nelemq,l)^(p-1);
                end
            end
            %Interior block
            for l=1:order*(order-1)+order-1
                lunk(3*(order+1)+l)=unkti(j,l);
                ldr(3*(order+1)+l)=1;
            end
        end
        etil=zeros(2,1);
        rcetil=0;
        icetil=0;
        for k=1:esize
            if lunk(k) ~=0
                etil(1)=etil(1)+(phxy(k,1)'.*(sol(lunk(k))).*ldr(k));
                etil(2)=etil(2)+(phxy(k,2)'.*(sol(lunk(k))).*ldr(k));
            end
        end
        
        out(splitt(j,q),1)=out(splitt(j,q),1)+etil(1);
        out(splitt(j,q),2)=out(splitt(j,q),2)+etil(2);
        help(splitt(j,q))=help(splitt(j,q))+1;
    end
end


for i=1:npoinn
    if help(i) ~= 0
        out(i,1)=out(i,1)./help(i);
        out(i,2)=out(i,2)./help(i);
    else
        disp('point not found')
    end
end


for i=1:2
    figure
    % call to plotting of values
    if i==1
        quiver(coordn(:,1),coordn(:,2),real(out(:,1)),real(out(:,2)))
        title('Quiver plot of Re(E^s)')
    else
        quiver(coordn(:,1),coordn(:,2),imag(out(:,1)),imag(out(:,2)))
        title('Quiver plot of Im(E^s)')
    end
end

