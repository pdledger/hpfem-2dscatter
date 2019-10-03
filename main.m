function main%Edge element solution of 2D scattering problemsclose alltic% add this directory and all subfoldersaddpath(genpath('./'))% define all data for this problem% problem = 1 PEC scattering by a cylinder% problem = 2 dielectric scattering by a cylinder% problem = 3 PEC Cylinder by a cylinder with grading factor in PML layer% problem = 4 dielectric scattering by a cylinder with grading factor in%             PML layer% problem = 5 PEC scattering by a NACA0012 problem=1;if problem==1probdata=problem1();elseif problem==2probdata=problem2();elseif problem==3probdata=problem3();elseif problem==4probdata=problem4();elseif problem==5probdata=problem5();% elseif problem==6% % probdata=problem6();% % elseif problem==7% % probdata=problem7();elsedisp('This problem file has not been created yet, press ctrl-c and start again')end% Set order of elements;order=probdata.order;% set mesh spacingh=probdata.h;% generate mesh using meshgenMesht=[];[nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals,bctype]=meshgen(Mesht,h,probdata);% Determine the number of edges[ne,edges,glob,dir]=edgeno(Mesht,Meshq,bsido,nelemt,nelemq,nboun,probdata);% Number the unknowns[order,unkl,unkh,unkqi,unkti,nunk,nunki,ndir,nunkpmlz,glob,edges,xmin,ymin,xmax,ymax,nedgepml]=nounk(ne,edges,nelemq,nelemt,glob,Mesht,Meshq,order,bctype,probdata);% enforce curved elements on the scatterer% default values -1 no curved edgeslingeom=probdata.lingeom;coefft=-1*ones(nelemt,3,2);coeff=-1*ones(nelemq,4,2);if lingeom==0     for i=1:length(probdata.rin)        rin=probdata.rin(i);        bcscatter=probdata.bcscatter(i);        [coefft,coeff,Mesht,Meshq]=getcoeff(nelemq,nelemt,Mesht,Meshq,edges,glob,ne,coefft,coeff,bcscatter,rin);     endendirn=[];icn=[];nz=0;map=[];% Assemble Quadrilateral elementsdisp('begin assembly');[rhs,basisx,basisy,curlbasis,x,w,rbasisx,rbasisy,dirval,tx,ty,I,J,X,Xabs,nz]=assembleq(nelemq,Meshq,...         order,map,irn,icn,nz,unkl,unkh,unkqi,...         unkti,glob,nunk,dir,nunki,edges,xmin,ymin,xmax,ymax,probdata,coefft,coeff);% Assemble Triangular Elementsdisp('begin assembly');[rhs,stiff,basisxt,basisyt,curlbast,intxi,inteta,intw,rbasisxt,rbasisyt,xt,wt,stiffabs,dirval]=assemblet(nelemt,Mesht,...         order,map,irn,icn,nz,unkl,unkh,...      	 unkti,glob,nunk,dir,nunki,edges,nelemq,x,w,rhs,dirval,xmin,ymin,tx,ty,probdata,coefft,coeff,I,J,X,Xabs);	 % Solve for Edge based Degree's of freedomdisp('Solve system');clear irn;clear icn;clear map;clear I;clear J;clear X;clear Xabs;sol=stiff\rhs;res=[];% introduce Dirichlet values back to solution vector[sol,unkl,unkh,nunki]=introd(sol,dirval,unkl,unkh,order,ne,nunki);disp('Assemble interiors (quads)'); [sol]=assembleqint(nelemq,Meshq,order,unkl,unkh,unkqi,...        unkti,glob,nunk,dir,nunki,edges,x,w,sol,basisx,basisy,curlbasis,...        rbasisx,rbasisy,xmin,ymin,tx,ty,probdata,coefft,coeff);disp('Assemble interiors (triangles)');if order >1[sol]=assembletint(nelemt,Mesht,...         order,unkl,unkh,unkti,glob,nunk,dir,nunki,edges,nelemq,xt,wt,...	 basisxt,basisyt,curlbast,intxi,inteta,intw,rbasisxt,rbasisyt,sol,...	 xmin,ymin,tx,ty,probdata,coefft,coeff);enddisp('Determine RCS');if nelemt==0[data]=rcs(sol,kx,ky,basisx,basisy,curlbasis,nelemq,Meshq,order,unkl,unkh,...         unkqi,unkti,glob,nunk,dir,nunki,x,w,Mesht,intxi,inteta,intw,...	     basisxt,basisyt,curlbast,nelemt,omega,edges,rbasisx,rbasisy,...	     probdata);else[data]=rcst(sol,basisx,basisy,curlbasis,nelemq,Meshq,order,unkl,unkh,...         unkqi,unkti,glob,nunk,dir,nunki,xt,wt,Mesht,intxi,inteta,intw,...	     basisxt,basisyt,curlbast,nelemt,edges,rbasisxt,rbasisyt,...	     probdata,coefft,coeff);	     end% Split the mesh for the quiver plot[Meshsq,Meshst,splitq,orders,splitt]=split(ne,0,Meshq,nelemq,Mesht,nelemt,dir,glob,edges,coefft,coeff);% Plot the quiver of the computed solutionmycontour(nelemt,nelemq,orders,Meshst,Meshsq,order,splitq,unkl,unkh...         ,unkqi,unkti,dir,glob,Meshq,Mesht,sol,splitt,coefft,coeff,probdata.omega);% Split the mesh[Meshsq,Meshst,splitq,orders,splitt]=split(ne,order,Meshq,nelemq,Mesht,nelemt,dir,glob,edges,coefft,coeff);% Plot the curl of the computed solutionmycontourz(nelemt,nelemq,orders,Meshst,Meshsq,order,splitq,unkl,unkh...         ,unkqi,unkti,dir,glob,Meshq,Mesht,sol,splitt,coefft,coeff,probdata.omega,xmin,ymin,xmax,ymax);     toc