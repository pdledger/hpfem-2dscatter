function varargout = plot_Mesh(Mesh)
% PLOT_MESH Mesh plot.
%
%   PLOT_MESH(MESH) generate 2D plot of the mesh.
%
%   H = PLOT_MESH(MESH) also returns the handle to the figure.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where M
%                is equal to the number of vertices contained in the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh. 
%
%   Example:
%
%   plot_Mesh(Mesh);
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constant
   
  OFFSET = 0.05;  % Offset parameter 
    
  % Compute axes limits
  
  X = Mesh.Coordinates(:,1);
  Y = Mesh.Coordinates(:,2);
  XMin = min(X);
  XMax = max(X);
  YMin = min(Y);
  YMax = max(Y);
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
    
  % Plot edges 
   npoin=length(X);
   U=zeros(1,npoin);
 % h = figure; 
  patch('faces',Mesh.Elements, ...
        'vertices',[X Y], ...
        'facecolor','none', ...
        'edgecolor','k');
  axis('equal');
  box('on');
  set(gca,'XLim',XLim,'YLim',YLim);
  title(['{\bf 2D mesh}']);
%  if(isfield(Mesh,'Edges'))
%    xlabel(['{\bf # Vertices  :  ', int2str(size(Mesh.Coordinates,1)), ...
%            ',      # Elements  :  ', int2str(size(Mesh.Elements,1)), ...
%            ',      # Edges  :  ',int2str(size(Mesh.Edges,1)),'}']);
%  elseif(isfield(Mesh,'BdEdges'))
%    xlabel(['{\bf # Vertices  :  ', int2str(size(Mesh.Coordinates,1)), ...
%            ',      # Elements  :  ', int2str(size(Mesh.Elements,1)), ...
%            ',      # Boundary edges  :  ',int2str(size(Mesh.BdEdges,1)),'}']);
%  else
%    xlabel(['{\bf # Vertices  :  ', int2str(size(Mesh.Coordinates,1)), ...
%            ',      # Elements  :  ', int2str(size(Mesh.Elements,1)),'}']); 
%  end
  drawnow;
   
  % Assign output arguments    
          
  if(nargout > 0)
    varargout{1} = h;
  end
  
return  


!DSPAM:42a01abf14392639321411!
