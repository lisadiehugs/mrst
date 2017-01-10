function plotCompareChooseBasis(basis, G, CG, Kfrac)
% plotCompareChooseBasis(basis, G, CG, Kfrac) plots the basis functions 
% contained in basis for different permeability values Kfrac
%
% REQUIRED PARAMETERS:
%
%   G  - Grid data structure
%
%   CG - coarse grid data structure
%
%   basis  -  cell structure containing one or several basis functions
%
%   ccb - index of coarse cell to compare basis function
%

% initially, plot basis for first coarse cell
ccb=1;
no_b = size(basis,1); % number of basis functions to plot
no_m = G.Matrix.cells.num; % number of cells in matrix grid
ifcells = CG.cells.interaction{ccb}; % fine cells in interaction region of coarse cell ccb
x = G.cells.centroids(1:no_m,1);
y = G.cells.centroids(1:no_m,2);
tri = delaunay(x,y);
hcb = zeros(no_b, 1); % checkbox handles
h = zeros(no_b, 1); % plot handles
figure(1); hold on;
% plot basis for coarse cell 1
for i=1:size(basis,1)
        z = full(basis{i}.B(1:no_m,ccb));
        h(i) = trimesh(tri,x,y,z);
end
axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1])
view(-20, 20);
colors = rand(no_b,3);
% set checkboxes to enable viewing of single basis functions
for i = 1:no_b
   hcb(i) = uicontrol('Style','checkbox','Value',1,...
                       'Position',[10 8+i*30 100 30],'String', ['K_f =' num2str(Kfrac(i))]);
   set(h(i), 'EdgeColor', 'k');
   set(h(i), 'FaceAlpha', 0.5); 
   set(h(i), 'FaceColor', colors(i,:))
end
set(hcb,'Callback',{@box_value,hcb,h});
% set popupmenu to choose coarse cell for which basis functions are plotted
strings = cell(CG.cells.num, 1);
for i = 1:CG.cells.num
    strings{i} = strcat(num2str(i));
end
hpop = uicontrol('Style', 'popup',...
           'String', strings,...
           'Position', [30 340 70 50]); 
set(hpop,'Callback',{@popup,h, basis, G, CG, tri});
end

% update figure (set new coarse cell)
function h = update(h, basis, G, CG, tri, ccb) 
no_m = G.Matrix.cells.num;
ifcells = CG.cells.interaction{ccb};
hold on;
    for i =1:size(basis,1)
        z = full(basis{i}.B(1:no_m,ccb));
        temp = trimesh(tri,G.cells.centroids(1:no_m,1),G.cells.centroids(1:no_m,2),z);
        set(h(i), 'ZData', temp.ZData);
        set(h(i), 'XData', temp.XData);
        set(h(i), 'YData', temp.YData);
        delete(temp);
    end
    refreshdata;
axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1])
end

%callback for popup menu
function popup(hObject, event, handles, basis, G, CG, tri)
index = get(hObject,'Value');
update(handles, basis, G, CG, tri, index);
end

%callback for checkboxes
function box_value(hObj,event,hcb,hax)
        % Called when boxes are used
        v = get(hObj,'Value');
        I = find(hcb==hObj);
        %[axes visibility]:
        s = {'off','on'};
        set(hax(I),'Visible',s{v+1})
end
