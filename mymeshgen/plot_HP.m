function varargout = plot_LFE(U,Mesh)
% PLOT_LFE Plot finite element solution.
%
%   PLOT_LFE(U,MESH) generates a plot fo the finite element solution U on the
%   mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where M is
%                equal to the number of vertices contained in the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M is
%                equal to the number of elements contained in the mesh.
%
%   H = PLOT(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_LFE(U,MESH);
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Generate figure
    
  %h = figure;
  patch('faces', Mesh.Elements, ...
        'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2)], ...
        'CData', U', ...
        'facecolor', 'interp', ...
        'edgecolor', 'k');
  axis('equal');
 
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = h;
  end
  
return

!DSPAM:42a01bbb14399383019716!
