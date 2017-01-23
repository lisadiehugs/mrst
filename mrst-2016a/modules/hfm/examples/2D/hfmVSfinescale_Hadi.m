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
celldim = [225 225];
physdim = [9 9];
celldim_hfm = [27 27];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
% define fracture lines 
fl = [4.5, 2 , 4.5, 7;
          2, 4.5, 7, 4.5];
hfmG = cartGrid(celldim_hfm, physdim);
hfmG = computeGeometry(hfmG);

% process fractures for hfm
[hfmG,fracture] = processFracture2D(hfmG, fl, 'verbose', mrstVerbose);
fracture.aperture = physdim(1)/celldim(1); % Fracture aperture

%% Compute CI and construct fracture grid
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
hfmG = CIcalculator2D(hfmG,fracture);
min_size = 0.3; cell_size = 0.32; % minimum and average cell size.
[hfmG,F,fracture] = gridFracture2D(hfmG,fracture,'min_size',min_size,'cell_size',cell_size);
figure; 
plotFractureNodes2D(hfmG,F,fracture);

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
fraccells = (G.cells.centroids(:,2)>(4.5-eps) & G.cells.centroids(:,2)<(4.5+eps) & G.cells.centroids(:,1)>2 & G.cells.centroids(:,1)<7 | ...
            G.cells.centroids(:,1)>(4.5-eps) & G.cells.centroids(:,1)<(4.5+eps) & G.cells.centroids(:,2)>2 & G.cells.centroids(:,2)<7);
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

%% Add BC

bc = [];
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = addBC(bc, left, 'pressure', 1*barsa, 'sat', [1 0]);
bc = addBC(bc, right, 'pressure', 0*barsa, 'sat', [0 1]);

bc_hfm = [];
xf = hfmG.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc_hfm = addBC(bc_hfm, left, 'pressure', 1*barsa, 'sat', [1 0]);
bc_hfm = addBC(bc_hfm, right, 'pressure', 0*barsa, 'sat', [0 1]);

%% Initialize state variables
% Once the wells are added, we can generate the components of the linear
% system corresponding to the two wells and initialize the solution
% structure (with correct bhp)
% fine-scale
state  = initResSol (G, 0*barsa, [0 1]);
%hfm
hfm_state  = initResSol (hfmG, 0*barsa, [0 1]);

%% Solve initial pressure in reservoir
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off
state = incompTPFA(state, G, hT, fluid, 'bc', bc);
hfm_state = incompTPFA(hfm_state, hfmG, T, fluid,  'bc', bc_hfm, 'use_trans',true);

%% Plot initial pressure

figure;
plotToolbar(G, state.pressure)
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
cx=caxis();
colormap jet(25)
axis tight off
title('Initial Pressure: Fine scale')

figure;
plotToolbar(hfmG, hfm_state.pressure)
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
caxis(cx);
view(90, 90)
axis tight off
title('Initial Pressure: HFM')
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
nt     = 30;
t200   = 0.1*(sum(pv)/sum(state.flux(left)));
Time   = t200;
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);
plotNo = 1;

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
figure;
while count <= nt,
   state = implicitTransport(state, G, dT, rock, fluid,  'bc',bc);
   hfm_state = implicitTransport(hfm_state, hfmG, dT, hfmG.rock, fluid,  'bc', bc_hfm, 'Trans', T);

   % Check for inconsistent saturations
   s = [state.s(:,1); hfm_state.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   state  = incompTPFA(state , G, hT, fluid,  'bc',bc);
   hfm_state  = incompTPFA(hfm_state, hfmG, T, fluid, 'bc', bc_hfm, 'use_trans',true);
   
   sol_fs{count,1} = state;
   sol_hfm{count,1} = hfm_state;
   
   %increase time
   t = t + dT;
   waitbar(t/Time,hwb);
   if ( t < plotNo*dTplot && t < Time), continue, end
    
    % Plot saturation
    heading = 'Saturation';
    
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.44]), cla
    plotCellData(G, state.s(:,1), 'edgealpha', 0.1);
    colormap(flipud(jet))
    view(0,90), axis equal off, title(['fine-scale t = ' num2str(round(t/day, 1)) ' days'] ,'FontSize',8)
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.44]), cla
    plotCellData(hfmG, hfm_state.s(:,1), 'edgealpha', 0.1);
    colormap(flipud(jet))
    view(0,90), axis equal off, title(['HFM t = ' num2str(round(t/day, 1)) ' days'],'FontSize',8)
    
    plotNo = plotNo+1;
    drawnow
   
   count = count + 1;
end
close(hwb);


    
%figure;
% subplot(1,2,1)
% plot(pvi,pfs(:,1),'-o', hfm_pvi, phfms, '-x');
% %plot(linspace(0,T,nt),pfs(:,1),'-o');
% leg = legend('Fine-scale', 'HFM', 'Location','Best');
% ylabel('Saturation at producer');
% xlabel('PVI [%]'); 
% set(gca,'XGrid','on','YGrid','on');
% axis tight
% subplot(1,2,2)
% plot(linspace(0,Time,nt),pfs(:,1),'-o', linspace(0,Time,nt), phfms, '-x');
% leg = legend('Fine-scale', 'HFM', 'Location','Best');
% ylabel('Saturation at producer');
% xlabel('Time [s]'); 
% set(gca,'XGrid','on','YGrid','on');
% axis tight
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
