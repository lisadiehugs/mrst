function varargout = plotCompareBasis(basis, ccb, G, CG, Kfrac)
% plotCompareBasis(basis, G) plots the basis functions contained in basis
%
% REQUIRED PARAMETERS:
%
%   G  - Grid data structure
%
%   p  -  cell structure containing one or several basis functions
%
%   ccb - index of coarse cell to compare basis function
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces.  OPTIONAL.

no_b = size(basis,1); % number of basis functions to plot
no_m = G.Matrix.cells.num; % number of cells in matrix grid
ifcells = CG.cells.interaction{ccb}; % fine cells in interaction region of coarse cell ccb
x = G.cells.centroids(1:no_m,1);
y = G.cells.centroids(1:no_m,2);
tri = delaunay(x,y);
hcb = zeros(no_b, 1); % checkbox handles
h = zeros(no_b, 1); % handles
figure; hold on;
for i =1:size(basis,1)
    z = full(basis{i}.B(1:no_m,ccb));
    h(i) = trimesh(tri,x,y,z); 
end
axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1])
view(-20, 20);
colors = rand(no_b,3);
for i = 1:no_b
   hcb(i) = uicontrol('Style','checkbox','Value',1,...
                       'Position',[10 8+i*30 100 30],'String', ['K_f =' num2str(Kfrac(i))]);
   set(h(i), 'EdgeColor', 'k');
   set(h(i), 'FaceAlpha', 0.5); 
   set(h(i), 'FaceColor', colors(i,:))
end
set(hcb,'Callback',{@box_value,hcb,h});
end

function box_value(hObj,event,hcb,hax)
        % Called when boxes are used
        v = get(hObj,'Value');
        I = find(hcb==hObj);
        %[axes visibility]:
        s = {'off','on'};
        set(hax(I),'Visible',s{v+1})
end