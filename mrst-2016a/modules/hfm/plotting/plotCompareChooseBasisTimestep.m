function plotCompareChooseBasisTimestep(basis, G, CG, paramval, paramstr)
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
nt = zeros(no_b);
for i=1:no_b
    nt(i) = size(basis{i},1);
end % number of time steps
no_m = G.Matrix.cells.num; % number of cells in matrix grid
ifcells = CG.cells.interaction{ccb}; % fine cells in interaction region of coarse cell ccb
x = G.cells.centroids(1:no_m,1);
y = G.cells.centroids(1:no_m,2);
tri = delaunay(x,y);
hcb = zeros(no_b, 1); % checkbox handles
h = zeros(no_b, 1); % plot handles
figure; hold on;
% plot basis for coarse cell 1 and first time step
for i=1:no_b
        z = full(basis{i}{1}.B(1:no_m,ccb));
        h(i) = trimesh(tri,x,y,z);
end
axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1])
view(-20, 20);
set(gcf,'Position',[200 200 1100 450]);
set(gca,'OuterPosition',[0.1 0 0.45 1]);
ax1 =  gca;
colors =          [0,        0.4470,  0.7410;
                  0.8500,    0.3250,   0.0980;
                  0.9290,    0.6940,   0.1250;
                  0.4940,    0.1840,   0.5560;
                  0.4660,    0.6740,   0.1880;
                  0.3010,    0.7450,   0.9330;
                  0.6350,    0.0780,   0.1840];
% set checkboxes to enable viewing of single basis functions
for i = 1:no_b
   hcb(i) = uicontrol('Style','checkbox','Value',1,...
                       'Position',[30 8+i*30 100 30],'String', [paramstr ' = ' num2str(paramval(i))]);
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
           'Position', [30 360 70 50]);
% Add a text uicontrol to label the slider.
txthpop = uicontrol('Style','text',...
            'Position',[30 410 70 20],...
            'String','CC Idx');

%set popupmenu to choose time step
stringsnt = cell(nt(1), 1);
for i = 1:nt(1)
    stringsnt{i} = strcat(num2str(i));
end
hpopnt = uicontrol('Style', 'popup',...
           'String', stringsnt,...
           'Position', [30 300 70 50]); 
txthpopnt = uicontrol('Style','text',...
            'Position',[30 350 70 20],...
            'String','Timestep');

        % plot grid and mark coarse cell
yy = subplot('position',[0.6 0.3 0.35 0.45]); hold on;
plotGrid(G, 'FaceColor', 'none', 'Edgecolor', [0.9,0.9,0.9]);
plotGrid(CG, 'FaceColor', 'none', 'Edgecolor', [0.4,0.4,0.4]);
axis tight;
hCCnode = zeros(CG.cells.num,1); 
hCCface = zeros(CG.cells.num,1);
hCCir= zeros(CG.cells.num,1);
for i = 1:CG.cells.num
    hCCir(i) = plotCellData(G, ones(size(CG.cells.interaction{i},1),1), CG.cells.interaction{i}, 'Visible', 'off', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.3);
    hCCface(i) = plotFaces(CG, CG.cells.faces(CG.cells.facePos(i):(CG.cells.facePos(i+1)-1)), 'r', 'Visible', 'off');
    hCCnode(i) = plot(CG.cells.centroids(i, 1),CG.cells.centroids(i, 2),'.r','markersize',12, 'Visible', 'off');
end
set(hCCir(1), 'Visible', 'on'); 
set(hCCnode(1), 'Visible', 'on');
set(hCCface(1), 'Visible', 'on');
% set callbacks
set(hpopnt,'Callback',{@popupnt,h, hpop, basis, G, CG, tri});
set(hpop,'Callback',{@popupccb,h, hpopnt, basis, G, CG, tri, ax1, hCCnode, hCCface, hCCir});
end

% update figure (set new coarse cell)
function h = updateccb(h, basis, G, CG, tri, ccb, nt, ax1, hCCnode, hCCface, hCCir) 
no_m = G.Matrix.cells.num;
ifcells = CG.cells.interaction{ccb};
hold on;
    for i =1:size(basis,1)
        z = full(basis{i}{nt}.B(1:no_m,ccb));
        temp = trimesh(tri,G.cells.centroids(1:no_m,1),G.cells.centroids(1:no_m,2),z);
        set(h(i), 'ZData', temp.ZData);
        set(h(i), 'XData', temp.XData);
        set(h(i), 'YData', temp.YData);
        delete(temp);
    end
    refreshdata;
    axis(ax1, [min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1])
    hCCnodeold= findobj(hCCnode, 'Visible', 'on');
    hCCfaceold= findobj(hCCface, 'Visible', 'on');
    hCCirold= findobj(hCCir, 'Visible', 'on');
    set(hCCnodeold, 'Visible', 'off');
    set(hCCfaceold, 'Visible', 'off');
    set(hCCirold, 'Visible', 'off');
    set(hCCnode(ccb), 'Visible', 'on');
    set(hCCface(ccb), 'Visible', 'on');
    set(hCCir(ccb), 'Visible', 'on');
end

% update figure (set new time step)
function h = updatetime(h, basis, G, CG, tri, nt, ccb) 
no_m = G.Matrix.cells.num;
ifcells = CG.cells.interaction{ccb};
hold on;
    for i =1:size(basis,1)
        z = full(basis{i}{nt}.B(1:no_m,ccb));
        temp = trimesh(tri,G.cells.centroids(1:no_m,1),G.cells.centroids(1:no_m,2),z);
        set(h(i), 'ZData', temp.ZData);
        set(h(i), 'XData', temp.XData);
        set(h(i), 'YData', temp.YData);
        delete(temp);
    end
    refreshdata;
%axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1])
end

%callback for popup menu for coarse cell
function popupccb(hObject, event, handles, hpopnt, basis, G, CG, tri, ax1, hCCnode, hCCface, hCCir)
ccb = get(hObject,'Value');
nt = get(hpopnt, 'Value');
updateccb(handles, basis, G, CG, tri, ccb, nt, ax1, hCCnode, hCCface, hCCir);
end

%callback for popup menu for time step
function popupnt(hObject, event, handles, hpop, basis, G, CG, tri)
nt = get(hObject,'Value');
ccb = get(hpop, 'Value');
updatetime(handles, basis, G, CG, tri, nt, ccb);
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
