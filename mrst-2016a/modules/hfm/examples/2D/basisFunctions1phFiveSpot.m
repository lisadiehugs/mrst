%% Introduction to the F-MsRSB Solver
% In this example, we will introduce you to the multiscale restriction
% smoothed basis (MsRSB) method for computing flow in embedded fracture
% models. We consider a 2D problem with no-flow boundary conditions and two 
% source terms at diagonally  opposite corners of a 2D Cartesian grid. 
%
% Notice that you need to have Metis installed to get this example to work.
% To get Metis working, you also need to set the global variable METISPATH.
% This can be done in your 'startup_user.m' file.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solversop
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path
    
%% Grid and fracture lines
% Construct a Cartesian grid comprising 90-by-90 cells, where each cell has
% dimension 1-by-1 m^2. Define 6 fracture lines.

%celldim = [99 99];
celldim = [18,18];
physdim = [99 99];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
% 3 horizontal and 3 vertical fractures (orthogonal) inside the domain
fl = [16.5, 0, 16.5, 99;
          49.5, 0, 49.5, 99;
          82.5, 0, 82.5, 99;
          0, 16.5, 99, 16.5;
          0, 49.5, 99, 49.5;
          0, 82.5, 99, 82.5];

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 2 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture
% figure(1);
% plotFractureLines(G, fracture);
% axis tight; 
% box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 0.05 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.8; cell_size = 1; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
% figure(1);
% plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix
% Set the permeability (K) as 1 Darcy in the matrix and 10000 Darcy in the
% fractures. Additionally, set the porosity of the matrix and fractures to
% 50%.
K_frac = [100;1000000];
basis_sb = cell(size(K_frac, 1), 1);
for i = 1:size(K_frac, 1)
dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy;
G.rock.poro = 0.5*ones(G.cells.num, 1);
poro_frac = 0.5;
G = makeRockFrac(G, K_frac(i), 'permtype','homogeneous','porosity',poro_frac);

%% Define fluid properties
% Define a single fluid of viscosity 1 cP and density 1000 kg/m3.

fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);

%% Define fracture connections as NNC and compute the transmissibilities
% In this section, we use the function defineNNCandTrans to combine the
% fracture and matrix grid structures into a single grid structure. In
% addition to that, we assign a 'non-neighbouring connection (NNC)' status
% to every fracture-matrix connection. For example, if fracture cell 'i' is
% located in the matrix cell 'j', then the connection i-j is defined as an
% NNC. To compute the flux between these elements, we compute a
% transmissibility for each NNC using the CI's computed earlier. Vector T
% contains the transmissibility for each face in the combined grid and each
% NNC.

[G,T] = defineNNCandTrans(G,F,fracture);

%% Add wells
% We will include two wells, both the injector and the producer controlled
% by bottom-hole pressure. The injector and producer are located at the
% bottom left and top right corner of the grid, respectively. Wells are
% described using a Peaceman model, giving an extra set of equations that
% need to be assembled, see simpleWellExample.m for more details on the
% Peaceman well model.

W = addWell([], G.Matrix, G.Matrix.rock, 1,'Type', ...
    'bhp', 'Val', 120*barsa,'Radius', .1);
W = addWell(W, G.Matrix, G.Matrix.rock, prod(celldim), 'Type', ...
    'bhp' , 'Val', 100*barsa, 'Radius', .1);

%% Initialize state variables
% Once the transmissibilities are computed, we can generate the
% transmissiblity matrix 'A' using the 'two-point flux approximation'
% scheme and initialize the solution structure.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);
% Get A matrix without source
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true);

%% Setup multiscale grids
% Next, we define a 9-by-9 matrix coarse grid and 4 coarse blocks (or
% coarse degrees of freedom) in the fracture. Additionally, we also define
% the support regions for the fracture and matrix basis functions. The
% matrix and fracture coarse grids are plotted.

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');

coarseDims = [6 6]; 
dof_frac = 12; %36; % Number of coarse blocks in the fracture grid
[CG, CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims,...
            'dof_frac',dof_frac);

%figure(1+(i-1)*3); plotFractureCoarseGrid2D(G,CG.partition,F);

%% Incompressible fine-scale solver
% The fine scale pressure solution is computed using the boundary
% conditions provided and the transmissiblity matrix computed earlier.

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
             'Wells',W, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions
% Using the transmissibility matrix 'A' we compute the basis functions for
% each fracture and matrix coarse block using the restriction smoothed
% basis method. Note that the matrix 'A' does not contain any source terms
% or boundary conditions. They are added to the coarse linear system when
% computing the multiscale pressure in the next section.

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb{i} = getMultiscaleBasis(CG, A, 'type', 'rsb');

%% Compute multiscale solution

% dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
% [state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb{i},...
%                'Wells', W,'use_trans',true);

%% Solve using MS-ILU and MS-GMRES
% Compute an iterative multiscale solution using ILU and GMRES
% preconditioners.

% fn = getSmootherFunction('type', 'ilu');
% 
% [~,report] = incompMultiscale(state, CG, T, fluid, basis_sb{i},...
%      'Wells', W, 'use_trans',true, 'tolerance', 1e-8, 'iterations', 100,...
%     'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);
% 
% [~,report2] = incompMultiscale(state, CG, T, fluid, basis_sb{i},...
%      'Wells', W, 'use_trans',true, 'tolerance', 1e-8, 'iterations', 100,...
%     'useGMRES', true, 'reconstruct', true, 'getSmoother', fn);

%% Plot results and convergence
% figure(2+(i-1)*3); colormap jet(25)
% plotCellData(G, state_fs.pressure,'EdgeColor','none')
% line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
% view(90, 90); colorbar, cx=caxis();
% axis tight off
% title(['Fine scale, K_f = ' num2str(K_frac(i))])
% 
% figure(3+(i-1)*3); colormap jet(25)
% plotCellData(G, state_ms.pressure,'EdgeColor','none')
% line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
% view(90, 90); colorbar, caxis(cx);
% axis tight off
% title(['F-MsRSB, K_f = ' num2str(K_frac(i))])

end

figure(size(K_frac, 1)*3+1);
plotToolbar(G,basis_sb{1}.B,'filterzero',true);
plotGrid(CG,'FaceColor','none');
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
axis tight; colorbar; 
title(['Basis Functions for fracture permeability: ' num2str(K_frac(1))]);

figure(size(K_frac, 1)*3+2);
plotToolbar(G,basis_sb{2}.B,'filterzero',true);
plotGrid(CG,'FaceColor','none');
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
axis tight; colorbar; 
title(['Basis Functions for fracture permeability: ' num2str(K_frac(2))])

%% Plot basis functions
Tri = delaunay(G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2));
ccb = 22; % index of coarse cell to compare basis function
ifcells = CG.cells.interaction{ccb};% fine cells in interaction region of coarse cell ccb
% plot basis functions in matrix cells for both fracture permeabilities and
% absolute deviation 
figure(size(K_frac, 1)*3+3);
pl1 = subplot(1,3,1);
%[xq,yq] = meshgrid(G.Matrix.cells.centroids(:,1), G.Matrix.cells.centroids(:,2));
%vq = griddata(G.Matrix.cells.centroids(:,1), G.Matrix.cells.centroids(:,2),full(basis_sb{1}.B(1:G.Matrix.cells.num,ccb)),xq,yq);
%surf(xq,yq,vq);
trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis_sb{1}.B(1:G.Matrix.cells.num,ccb)));
outlineCoarseCell(G, CG.partition, ccb, 'b');
title(['K_f = ' num2str(K_frac(1))])
colormap(pl1,spring)
pl2 = subplot(1,3,2);
trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis_sb{2}.B(1:G.Matrix.cells.num,ccb)));
colormap(pl2,winter)
title(['K_f = ' num2str(K_frac(2))])
axis([pl1, pl2],[min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1]); 
pl3 = subplot(1,3,3);
trimesh(Tri, G.cells.centroids(1:G.Matrix.cells.num,1), G.cells.centroids(1:G.Matrix.cells.num,2), full(basis_sb{1}.B(1:G.Matrix.cells.num,ccb))-full(basis_sb{2}.B(1:G.Matrix.cells.num,ccb)));
axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) -1 1]); 

% % plot basis functions in fracture cells for both fracture permeabilities and
% % absolute deviation 
% fstart = G.Matrix.cells.num+1; %first fracture cell in assembled grid
% Tri = delaunay(G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2));
% figure(size(K_frac, 1)*3+4);
% pl1 = subplot(1,3,1);
% %trimesh(Tri, G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2), full(basis_sb{1}.B(fstart:end,ccb)));
% plot3(G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2), full(basis_sb{1}.B(fstart:end,ccb)));
% title(['K_f = ' num2str(K_frac(1))])
% colormap(pl1,spring)
% pl2 = subplot(1,3,2);
% %trimesh(Tri, G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2), full(basis_sb{2}.B(fstart:end,ccb)));
% plot3(G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2), full(basis_sb{2}.B(fstart:end,ccb)));
% colormap(pl2,winter)
% title(['K_f = ' num2str(K_frac(2))])
% axis([pl1, pl2],[min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) 0 1]); 
% pl3 = subplot(1,3,3);
% %trimesh(Tri, G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2), full(basis_sb{1}.B(fstart:end,ccb))-full(basis_sb{2}.B(fstart:end,ccb)));
% plot3(G.cells.centroids(fstart:end,1), G.cells.centroids(fstart:end,2), full(basis_sb{1}.B(fstart:end,ccb))-full(basis_sb{2}.B(fstart:end,ccb)));
%axis([min(G.cells.centroids(ifcells,1)) max(G.cells.centroids(ifcells,1)) min(G.cells.centroids(ifcells,2)) max(G.cells.centroids(ifcells,2)) -1 1]); 


%% calculate deviation of basis functions in L2 
deviation = abs(full(basis_sb{1}.B(ifcells,ccb))-full(basis_sb{2}.B(ifcells,ccb)));
l2deviation = sqrt(sum(deviation.*deviation.*G.cells.volumes(ifcells),1)/sum(G.cells.volumes(ifcells),1));
l1deviation = sum(deviation.*G.cells.volumes(ifcells),1)/sum(G.cells.volumes(ifcells),1);


%%% Plot convergence
% figure;
% semilogy(report.resvec, ':+'); hold on;
% semilogy(report2.resvec, 's-');
% legend('Convergence of ILU(0)', 'Convergence of GMRES');
% title('Iterative convergence');

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