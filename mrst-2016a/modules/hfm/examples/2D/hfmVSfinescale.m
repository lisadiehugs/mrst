%% Two-Phase Problem with a Quarter Five-Spot Well Pattern
% test case to compare hfm with fine-scale resolution of fractures
mrstModule add incomp
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Define geometry and rock parameters
% Construct a Cartesian grid of size 297-297 cells. Set the permeability K 
% to be homogeneous, isotropic and equal 1 Darcy in the matrix and 10e4 Darcy 
% in the fractures. The porosity is 0.2 in the fractures and 0.5 in the matrix. 
% Compute one-sided (half) transmissibilities from input grid and rock properties.
celldim = [297 297];
physdim = [99 99];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
% define fracture lines 
fl = [16.5, 0, 16.5, 99;
          49.5, 0, 49.5, 99;
          82.5, 0, 82.5, 99;
          0, 16.5, 99, 16.5;
          0, 49.5, 99, 49.5;
          0, 82.5, 99, 82.5];
hfmG = cartGrid(celldim, physdim);
hfmG = computeGeometry(hfmG);

% process fractures for hfm
[hfmG,fracture] = processFracture2D(hfmG, fl, 'verbose', mrstVerbose);
fracture.aperture = physdim(1)/celldim(1); % Fracture aperture

%% Compute CI and construct fracture grid
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
hfmG = CIcalculator2D(hfmG,fracture);
min_size = 0.4; cell_size = 0.5; % minimum and average cell size.
[hfmG,F,fracture] = gridFracture2D(hfmG,fracture,'min_size',min_size,'cell_size',cell_size);

%% Set rock properties
K_frac = 10000; K_matrix = 1; % * Darcy
p_frac = 0.2; p_matrix = 0.5;
% hfm
hfmG.rock.poro = p_matrix*ones(hfmG.cells.num, 1);
hfmG.rock.perm = K_matrix*ones(hfmG.cells.num, 1)*darcy;
hfmG = makeRockFrac(hfmG, K_frac, 'permtype','homogeneous','porosity',p_frac);
% fine-scale
rock.perm  = ones(G.cells.num,1)*K_matrix *darcy;
rock.poro  = ones(G.cells.num, 1)*p_matrix;
eps = 1e-3;
fraccells = (G.cells.centroids(:,2)>(16.5-eps) & G.cells.centroids(:,2)<(16.5+eps) | ...
            G.cells.centroids(:,1)>(16.5-eps) & G.cells.centroids(:,1)<(16.5+eps) | ...
            G.cells.centroids(:,2)>(49.5-eps) & G.cells.centroids(:,2)<(49.5+eps) | ...
            G.cells.centroids(:,1)>(49.5-eps) & G.cells.centroids(:,1)<(49.5+eps) | ...
            G.cells.centroids(:,2)>(82.5-eps) & G.cells.centroids(:,2)<(82.5+eps) | ...
            G.cells.centroids(:,1)>(82.5-eps) & G.cells.centroids(:,1)<(82.5+eps));
rock.perm(fraccells) = K_frac * darcy;
rock.poro(fraccells) = p_frac;
rock.ntg   = ones(G.cells.num, 1);

%compute transmissibilities
%fine-scale
hT  = computeTrans(G, rock, 'Verbose', true);
figure;
plotCellData(G, rock.perm);
%hfm
[hfmG,T] = defineNNCandTrans(hfmG,F,fracture);

%% Define the two-phase fluid model
% The <matlab:help('initSimpleFluid') two-phase fluid model> has default values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 10] cP.
    fluid = initSimpleFluid('mu' , [   1,  2] .* centi*poise     , ...
            'rho', [1000, 700] .* kilogram/meter^3, ...
            'n'  , [   2,   2]);

%% Introduce wells
% We will include two wells, one rate-controlled vertical well and one
% horizontal well controlled by bottom-hole pressure. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled, see the <matlab:edit('incompTutorialWells.m') tutorial on well
% models> for more details.
W = addWell([], G, rock, 1,'InnerProduct', 'ip_tpf','Type', ...
    'bhp', 'Val', 120*barsa,'Radius', .06, 'Comp_i', [1, 0]);
    W = addWell(W, G, rock, prod(celldim), 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 100*barsa, 'Radius', .06, 'Comp_i', [0, 1]);

hfmW = addWell([], hfmG.Matrix, hfmG.Matrix.rock, 1,'InnerProduct', 'ip_tpf','Type', ...
    'bhp', 'Val', 120*barsa,'Radius', .06, 'Comp_i', [1, 0]);
    hfmW = addWell(hfmW, hfmG.Matrix, hfmG.Matrix.rock, prod(physdim), 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 100*barsa, 'Radius', .06, 'Comp_i', [0, 1]);

%% Initialize state variables
% Once the wells are added, we can generate the components of the linear
% system corresponding to the two wells and initialize the solution
% structure (with correct bhp)
% fine-scale
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);
%hfm
hfm_state  = initResSol (hfmG, 0);
hfm_state.wellSol = initWellSol(hfmW, 0);

%% Solve initial pressure in reservoir
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off
state = incompTPFA(state, G, hT, fluid, 'wells', W);
hfm_state = incompTPFA(hfm_state, hfmG, T, fluid, 'wells', hfmW, 'use_trans',true);

%% Transport loop
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 15 equally spaced time steps). The error introduced by this
% splitting of flow and transport can be reduced by iterating each time
% step until e.g., the residual is below a certain user-prescribed
% threshold (this is not done herein).
pv     = poreVolume(G,rock);
hfm_pv = poreVolume(hfmG,hfmG.rock);
nt     = 45;
Time      = 0.5*(sum(pv)/state.wellSol(1).flux);
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

%%
% The transport equation will be solved by the single-point upstream method
% with either explicit or implicit time discretizations. Both schemes may
% use internal time steps to obtain a stable discretization. To represent
% the two solutions, we create new solution objects to be used by the
% solver with implicit transport step.

%% Start the main loop
pvi = zeros(nt,1); hfm_pvi = zeros(nt,1);
pfs = zeros(nt,1); phfms = zeros(nt,1); 
volflux = zeros(nt,1); hfm_volflux = zeros(nt,1);
t  = 0;
hwb = waitbar(0,'Time loop');
count = 1;
sol_fs = cell(nt,1);
sol_hfm = cell(nt,1);
while t < Time,
   state = implicitTransport(state, G, dT, rock, fluid, 'wells', W);
   hfm_state = implicitTransport(hfm_state, hfmG, dT, hfmG.rock, fluid, 'wells', hfmW, 'Trans', T);

   % Check for inconsistent saturations
   s = [state.s(:,1); hfm_state.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   state  = incompTPFA(state , G, hT, fluid, 'wells', W);
   hfm_state  = incompTPFA(hfm_state, hfmG, T, fluid, 'wells', hfmW, 'use_trans',true);
   
   sol_fs{count,1} = state;
   sol_hfm{count,1} = hfm_state;
   
   %increase time
   t = t + dT;
   waitbar(t/Time,hwb);
   % pore volume injected
   volflux(count) = state.wellSol(1).flux*dT;
   hfm_volflux(count) = hfm_state.wellSol(1).flux*dT;
   pvi(count) = 100*(sum(volflux(1:count)))/sum(pv);
   hfm_pvi(count) = 100*(sum(hfm_volflux(1:count)))/sum(pv);
   pfs(count,1) = state.s(W(2).cells,1); 
   phfms(count,1) = hfm_state.s(hfmW(2).cells,1); 
   count = count + 1;
end
close(hwb);
%% Plot saturations
plotNo = 1;
figure; hold on; colormap(flipud(winter))
boundary = any(G.faces.neighbors==0,2);
hfmboundary = any(hfmG.faces.neighbors==0,2);
facelist = 1:G.faces.num;
hfmfacelist = 1:hfmG.faces.num;
bfaces = facelist(boundary);
hfmbfaces = hfmfacelist(hfmboundary);
for i = nt/3:nt/3:nt
    state_fs = sol_fs{i,1}; state_hfm = sol_hfm{i,1};
    heading = [num2str(round(pvi(i))),  ' % PVI'];
    % Plot saturation
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
    plotCellData(G,state_fs.s(:,1),'EdgeColor','none');
    line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
    plotFaces(G,bfaces,'k','linewidth',1)
    axis square off, 
    title(['Reference: ', heading],'FontSize',8);
    view(0,90); caxis([0 1]);
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
    plotCellData(hfmG,state_hfm.s(:,1),'EdgeColor','none');
    line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
    plotFaces(hfmG,hfmbfaces,'k','linewidth',1)
    axis square off, 
    title(['HFM: ',  heading],'FontSize',8);
    view(0,90); caxis([0 1]);
    plotNo = plotNo+1;
end
    
figure;
plot(pvi,pfs(:,1),'-o', hfm_pvi, phfms, '-x');
%plot(linspace(0,T,nt),pfs(:,1),'-o');
leg = legend('Fine-scale', 'HFM', 'Location','Best');
ylabel('Saturation at producer');
xlabel('PVI [%]'); 
set(gca,'XGrid','on','YGrid','on');
axis tight
%%

% <html>
% <p><font size="-1">
% Copyright 2009-2016 SINTEF ICT, Applied Mathematics.
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
