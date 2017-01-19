%% F-MsRSB applied to a heterogeneous 2D domain
% Two-phase example modeling water injection through a quarter 5-spot into a
% 2-dimensional fractured porous media. The flow problem is solved both by a
% fine-scale and a multiscale solver.

% Notice that you need to have Metis installed to get this example to work.
% To get Metis working, you also need to set the global variable METISPATH.
% This can be done in your 'startup_user.m' file.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Call run_compareParams_2phFiveSpot
K_frac = 1e2;
%paramstr = 'tolerance'; %'fracSupportRadius'
paramstr = 'fracSupportRadius';
%paramval = [5e-4; 5e-3; 5e-2; 5e-1];
paramval = [1;2;4;8];
basis = cell(size(paramval, 1), 1);
pv_out = cell(size(paramval, 1), 1); % pvi, pfs, pms
for i = 1:size(paramval, 1)
    [G, CG, basis{i}, pv_out{i}] = run_compareParams_2phFiveSpot(K_frac, paramstr,  paramval, i);
    
    % Plot water saturation at producer 
    figure;
    plot(pv_out{i}(:,1) ,pv_out{i}(:,2),'-o',pv_out{i}(:,1),pv_out{i}(:,3),'--*');
    leg = legend('Fine-scale','Multiscale','Location','Best');
    ylabel('Saturation at producer');
    xlabel('PVI [%]'); 
    set(gca,'XGrid','on','YGrid','on');
    axis([0 max(pv_out{i}(:,1)) 0 1]);
end
 
%% Plot basis functions
Tri = delaunay(G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2));
ccb = 14; % index of coarse cell to compare basis function
ifcells = CG.cells.interaction{ccb};% fine cells in interaction region of coarse cell ccb
% plot basis functions in matrix cells for both fracture permeabilities and absolute deviation 
% if (size(paramval, 1) == 1)
%     figure;
%     trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis{1}.B(1:G.Matrix.cells.num,ccb)));
%     outlineCoarseCell(G, CG.partition, ccb, 'b');
%     title(['K_f = ' num2str(K_frac(1))])
%     axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) -1 1]);
% elseif (size(paramval, 1) == 2)
%     figure;
%     pl1 = subplot(1,3,1);
%     trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis{1}.B(1:G.Matrix.cells.num,ccb)));
%     outlineCoarseCell(G, CG.partition, ccb, 'b');
%     title(['K_f = ' num2str(K_frac(1))])
%     colormap(pl1,spring(10))
%     pl2 = subplot(1,3,2);
%     trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis{2}.B(1:G.Matrix.cells.num,ccb)));
%     colormap(pl2,winter(10))
%     title(['K_f = ' num2str(K_frac(2))])
%     axis([pl1, pl2],[min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1]); 
%     pl3 = subplot(1,3,3);
%     trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis{1}.B(1:G.Matrix.cells.num,ccb))-full(basis{2}.B(1:G.Matrix.cells.num,ccb)));
%     axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) -1 1]); 
% end

plotCompareChooseBasis(basis, G, CG, paramval, paramstr);

% figure; hold on;
% for i = 1:size(K_frac, 1)
%     c(1:G.Matrix.cells.num) = i;
%     trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis{i}.B(1:G.Matrix.cells.num,ccb)), c);
%     axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1]); 
% end
% <html>
% <p><font size="-1">
% Copyright 2009-2016 TU Delft and SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>