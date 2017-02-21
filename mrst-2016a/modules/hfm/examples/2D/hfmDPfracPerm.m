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

%% Grid and fracture lines
% Construct a Cartesian grid comprising 50-by-20 cells, where each cell has
% dimension 10-by-10 m^2. Define 3 fracture lines by supplying their end
% points in the form [x1 y1 x2 y2].

 celldim = [120 70];
physdim = [100 50];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [0, 25, 25, 0;
      0, 25, 25, 50;
      75, 0, 100, 25;
      75, 50, 100, 25;

      % "inner" fractures 
      0,0,50,50;
      25, 0, 75, 50;
      50, 0, 100, 50;
      0, 50, 50, 0;
      25, 50, 75, 0;
      50, 50, 100, 0];
%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 3 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G, fl, 'verbose', mrstVerbose);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture);
axis equal tight; 
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 10 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.8; cell_size = 1; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); 
axis equal tight; box on

%% Set rock properties in fracture and matrix
% For this example, we will generate the porosity as a Gaussian field. To
% get a crude approximation to the permeability-porosity relationship, we
% assume that our medium is made up of uniform spherical grains of diameter
% dp = 10 m, for which the specic surface area is Av = 6 = dp. With these
% assumptions, using the Carman Kozeny relation, we can then calculate the
% isotropic permeability (K). The rock properties are then plotted.

K_frac = 1e4;
K_m = 1e-10;
dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = K_m*ones(G.cells.num,1)*darcy;
G.rock.poro = 0.4*ones(G.cells.num, 1); %0.2 works also fine
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
        'n'  , [   1,   1]);

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

%% Add BC
% Set boundary condition P = 10 bars on the left face and P = 1 bar on the
% right face of the domain.

bc = [];
xf = G.Matrix.faces.centroids(:, 1);
yf = G.Matrix.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = addBC(bc, left, 'pressure', 10*barsa, 'sat', [1 0]);
bc = addBC(bc, right, 'pressure', 1*barsa, 'sat', [0 1]);

boundfaces=findfracboundaryfaces2d(G,1e-5);
bc = addBC(bc,boundfaces.West,'pressure',10*barsa, 'sat', [1 0]);
bc = addBC(bc,boundfaces.East,'pressure',1*barsa, 'sat', [0 1]);

%% define outflow faces
    OutflowfacesMat = right;
    OutflowfacesFrac = boundfaces.East';
    Outflowfaces = [OutflowfacesMat;OutflowfacesFrac];

%% Initialize state variables
% Here, we initialize the solution structure using the combined grid and
% the wells defined above.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 1*barsa, [0 1]);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state = incompTPFA(state, G, T, fluid, 'bc', bc, 'MatrixOutput', true, 'use_trans',true);

%% Incompressible Two-Phase Flow
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 60 equally spaced time steps amounting to 50 % PV Injected).
% The error introduced by this splitting of flow and transport can be
% reduced by iterating each time step until e.g., the residual is below a
% certain user-prescribed threshold (this is not done herein).

nt     = 50; 
Time   = 5*day; 
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

sout = zeros(nt,1); % stores accumulated saturation along outflow boundary
smat = cell(nt,1); % stores accumulated saturation for domain
sfrac = cell(nt,1);
totflux = zeros(nt,1);
pvi = zeros(nt,1);
sol = cell(nt,1);
plotNo = 1;

t  = 0;
count = 1; 
while t < Time
    state = implicitTransport(state, G, dT, G.rock, fluid, 'bc', bc, 'Trans', T,'verbose',false);
    
    % Check for inconsistent saturations
    assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);
 
    
    % Update solution of pressure equation.
    state  = incompTPFA(state, G, T, fluid, 'bc', bc, 'use_trans',true);
    sol{count,1} = state;
    
    % sum up outflow
    totflux(count,1) = sum(state.flux(Outflowfaces),1);
    % saturation in domain for matrix (accumulate values into nmx blocks)
    nmx = 8;
    em = 0:physdim(1)/nmx:physdim(1);
    ym = discretize(G.cells.centroids(1:G.Matrix.cells.num, 1),em);
    smat{count} = accumarray(ym, state.s(1:G.Matrix.cells.num,1))./histcounts(G.cells.centroids(1:G.Matrix.cells.num, 1),em)';
    % saturation in domain for fracture (accumulate values into nfx blocks)
    nfx = 20;
    ef = 0:physdim(1)/nfx:physdim(1);
    yf=discretize(G.cells.centroids((G.Matrix.cells.num+1):end, 1), ef);
    sfrac{count} = accumarray(yf, state.s(G.Matrix.cells.num+1:end,1))./ histcounts(G.cells.centroids(G.Matrix.cells.num+1:end, 1),ef)';
    % integrate and store saturation at outflow boundary
    sout(count) = smat{count}(end);
    
    % plot saturation in domain
    if (mod(count,2)~=0)
        figure(4)
        set(4, 'Position', [100, 100, 900, 400])
        subplot(1,2,1)
        plot(em(1:end-1)+em(2)/2, smat{count}, '-o');
        ylabel('Saturation [-]');
        xlabel('x [m]');
        legend('Fine-scale','Location','Best');
        axis([0 physdim(1) 0 1])
        title(['Saturation in Matrix (' num2str(nmx) ' cells) - t = ' num2str(t/day) ' days'])
        subplot(1,2,2)
        plot(ef(1:end-1)+ef(2)/2, sfrac{count}, '-o');
        legend('Fine-scale','Location','Best');
        ylabel('Saturation [-]');
        xlabel('x [m]');
        axis([0 physdim(1) 0 1])
        title(['Saturation in Fractures (' num2str(nfx) ' cells) - t = ' num2str(t/day) ' days'])
    end
    %% increase time
    t = t + dT;
    count = count + 1;
    
    if ( t < plotNo*dTplot && t < Time), continue, end
        figure(5);
        set(5, 'Position', [100, 100, 1100, 500])
        % Plot saturation
        heading = [num2str(convertTo(t,day)),  ' days'];

        r = 0.01;
        subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.44]), cla
        plotCellData(G, state.s(:,1), 'edgealpha', 0.1);
        colormap(flipud(jet))
        view(0,90), axis equal off, title('saturation','FontSize',8)


        plotNo = plotNo+1;
        drawnow
    
end

 %% plot saturation and flux at outflow
    figure(7);
    subplot(1,2,1)
    plot(linspace(0, Time, nt)/day, sout,'-o');
    leg = legend('Fine-scale','Location','Best');
    ylabel('Saturation at Outflow boundary [-]');
    xlabel('Time [days]');
    %saveas(7, ['watercutKf=' num2str(K_frac) 'Km=' num2str(K_m) '.epsc']);
    subplot(1,2,2);
    plot(linspace(0, Time, nt)/day, totflux,'-o');
    leg = legend('Fine-scale','Location','Best');
    ylabel('Flux at Outflow boundary [-]');
    xlabel('Time [days]');
    

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
