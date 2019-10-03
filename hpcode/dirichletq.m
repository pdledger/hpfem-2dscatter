function [dirval]=dirichletq(edges,x,w,glob,coord,intmaq,unkl,unkh,dir,...
                             rbasisx,rbasisy,nelemq,order,coeff)

nip=length(x);

esize=4*(order+1)+2*order*(order+1);

norml=[0 -1; 0 1; 1 0; -1 0]';
tanl=[1 0; -1 0; 0 1; 0 -1]';
dirval=[];
for i=1:nelemq

    for l=1:4
        for k=1:2
            xyq(l,k)=coord(intmaq(i,l),k);
        end
    end

   for l=1:4
      for k=1:2
      localblend(l,k)=coeff(i,l,k);
      end
   end


    for jj=1:4
        if edges(glob(i,jj),3)~=0 & bctype(edges(glob(i,jj),3))==2
            %  then we have a Dirichlet edge
            num=0;
            a=zeros(order+1,order+1);
            r=zeros(order+1,1);
            den=zeros(order+1,1);
            temp=zeros(esize,1);

            for p=1:nip
                if jj==1
                    xi=x(p);
                    et=0;
                elseif jj==2
                    xi=x(p);
                    et=1;
                elseif jj==3
                    xi=1;
                    et=x(p);
                else
                    xi=0;
                    et=x(p);
                end
                phx=rbasisx(:,p+(nip*(jj-1)));
                phy=rbasisy(:,p+(nip*(jj-1)));
                ph=[phx phy];


                J=mapquad(xyq,xi,et,localblend);
                detJ=det(J);
                Jinv=(inv(J))';
                phxy=Jinv*ph';
                phxy=phxy';

                nm=inv(J)'*norml(:,jj);
                nm=nm./norm(nm);
                tang=J*tanl(:,jj);
                tang=tang./norm(tang);

                if jj==1 || jj==2
                    v=sqrt(J(1,1)^2+J(2,1)^2);
                else
                    v=sqrt(J(1,2)^2+J(2,2)^2);
                end
                [xp,yp]=getcoord(xyq,xi,et,localblend);


                %  use the problem file
                fun=probdata.dirfun;
                arg=probdata.dirfunarg;
                index=edges(glob(i,jj),3);
                tange=fun(xp,yp,index,nm,arg);

                %  For scattering problems
                %   e=[-sin(theta); cos(theta)]*(cos(omega*((cos(theta)*xp)+(sin(theta)*yp)))-...
                %   j*sin(omega*((cos(theta)*xp)+(sin(theta)*yp))));;

                %  the tangential components of all basis functions except those
                %  on edge j vanish
                for q=1:order+1
                    %   For Wave Propagation Problems
                    %    r(q)=r(q)+w(p)*v*(tang(1)*e(1)+tang(2)*e(2))*(nm(1)*phxy(jj+(4*(q-1)),2)-nm(2)*phxy(jj+(4*(q-1)),1));;
                    %   For Scattering Problems n x E= - nx E^in
                    r(q)=r(q)+w(p)*v*tange*(nm(1)*phxy(jj+(4*(q-1)),2)-nm(2)*phxy(jj+(4*(q-1)),1));;
                end
                for q=1:order+1
                    for qq=1:order+1
                        a(q,qq)=a(q,qq)+(w(p)*v*(nm(1)*phxy(jj+(4*(q-1)),2)-nm(2)*phxy(jj+(4*(q-1)),1))*...
                            (nm(1)*phxy(jj+(4*(qq-1)),2)-nm(2)*phxy(jj+(4*(qq-1)),1)));
                    end
                end

            end
            if det(a)==0
                disp('error we cannot apply Dirichlet BCs')
            end
            s=a\r;
            if dir(i,jj) < 0
                disp('wrong orientation defined in mesh!')
            end
            %  lowest order block
            dirval(abs(unkl(glob(i,jj))))=s(1)*dir(i,jj);
            %  higher order block
            for q=1:order
                dirval(abs(unkh(glob(i,jj),q)))=s(q+1)*(dir(i,jj)^(q-1));
            end
        end
    end
end
