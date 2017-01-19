%% Two-Phase Problem with a Quarter Five-Spot Well Pattern
% In this second introductory example to the HFM module, we show the impact
% of fractures on fluid migration using the hierarchical/embedded fracture
% model. To this end, we consider a two-phase example with three
% intersecting fractures in the center of the model. Oil is recovered by a
% production well in the NE corner, which is supported by a water-injector
% in the SW corner. 

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

% measure cputime
tic;

%% Grid and fracture line(s)
% Construct a Cartesian grid comprising 100-by-100 cells, where each cell
% has dimension 1-by-1 m^2. Define 1 fracture line, roughly 85 m in length,
% extending diagonally through the centre of the domain.

celldim = [99 99];
physdim = [99 99];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [16.5, 0, 16.5, 99;
          49.5, 0, 49.5, 99;
          82.5, 0, 82.5, 99;
          0, 16.5, 99, 16.5;
          0, 49.5, 99, 49.5;
          0, 82.5, 99, 82.5];

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 3 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

%dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G, fl, 'verbose', mrstVerbose);
fracture.aperture = 1/25; % Fracture aperture
% figure;
% plotFractureLines(G,fracture);
% axis equal tight; 
% box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 10 m. Next, the fracture grid is plotted on top of the matrix grid.

%dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.8; cell_size = 1; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
% clf; plotFractureNodes2D(G,F,fracture); 
% axis equal tight; box on

%% Set rock properties in fracture and matrix
% For this example, we will generate the porosity as a Gaussian field. To
% get a crude approximation to the permeability-porosity relationship, we
% assume that our medium is made up of uniform spherical grains of diameter
% dp = 10 m, for which the specic surface area is Av = 6 = dp. With these
% assumptions, using the Carman Kozeny relation, we can then calculate the
% isotropic permeability (K). The rock properties are then plotted.

%dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy;
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 1e2; % * Darcy in makeRockFrac
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','porosity',poro_frac);

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 1] cP.
% * corey-coefficient: [2, 2] = [2 2].

fluid = initSimpleFluid('mu' , [   1,  2] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

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
% We will include two wells, one rate-controlled injector well and the
% producer controlled by bottom-hole pressure. The injector and producer
% are located at the bottom left and top right corner of the grid,
% respectively. Wells are described using a Peaceman model, giving an extra
% set of equations that need to be assembled, see simpleWellExample.m for
% more details on the Peaceman well model.

W = addWell([], G.Matrix, G.Matrix.rock, 1,'InnerProduct', 'ip_tpf','Type', ...
    'bhp', 'Val', 120*barsa,'Radius', .1, 'Comp_i', [1, 0]);
W = addWell(W, G.Matrix, G.Matrix.rock, prod(celldim), 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 100*barsa, 'Radius', .1, 'Comp_i', [0, 1]);

%% Initialize state variables
% Here, we initialize the solution structure using the combined grid and
% the wells defined above.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state = incompTPFA(state, G, T, fluid, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);

% %% Plot initial pressure
% clf
% hp = plotCellData(G, state.pressure,'EdgeColor','none');
% xw = G.cells.centroids(vertcat(W.cells),:);
% hold on
% plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
% plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
% hold off
% colormap jet
% view(0, 90); colorbar; axis equal tight;
% title('Initial pressure');

%% Incompressible Two-Phase Flow
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 60 equally spaced time steps amounting to 50 % PV Injected).
% The error introduced by this splitting of flow and transport can be
% reduced by iterating each time step until e.g., the residual is below a
% certain user-prescribed threshold (this is not done herein).

pv     = poreVolume(G,G.rock);
nt     = 45;
Time   = 0.5*(sum(pv)/state.wellSol(1).flux);
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

pvi = zeros(nt,1);
sol = cell(nt,1);

t  = 0;
count = 1; title('Saturation');
%colorbar off; colormap(flipud(winter));
while t < Time,
    state = implicitTransport(state, G, dT, G.rock, fluid, 'wells', W, 'Trans', T,'verbose',true);
    
    % Check for inconsistent saturations
    assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);

    % Plot saturation
    %delete(hp)
    %hp = plotCellData(G,state.s,state.s>0); drawnow; pause(.1);
    
    % Update solution of pressure equation.
    state  = incompTPFA(state, G, T, fluid, 'wells', W, 'use_trans',true);
    sol{count,1} = state;
    
    t = t + dT;
    pvi(count) = 100*(sum(state.wellSol(1).flux)*t)/sum(pv);
    count = count + 1;
    
end

% %% Plot saturations
% clf, plotToolbar(G,cellfun(@(x) struct('s', x.s(:,1)), sol,'UniformOutput', false));
% line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
% hold on
% plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
% plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
% hold off
% colormap(flipud(winter)); caxis([0 1]); axis equal tight;
toc;
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
