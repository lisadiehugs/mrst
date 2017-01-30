%% quarterSpotDP - This example outlines the basic usage of the 
% two-phase dual-porosity model. We inject the first phase on 
% one corner and produce at the opposite one.
clc;
close all; clear all;

%% Include the essential modules. dual-porosity model has the functionality
% we need.
mrstModule add ad-props ad-core ad-blackoil ad-fi blackoil-sequential dual-porosity

%% Create the grid: a simple cartesian grid, 100m x 100m and 50 blocks on x
% and y directions.

% Two layer grid
y_size = 9;
x_size = 27;

Ny = 18;
Nx = 54;


G = cartGrid([Nx Ny],[x_size y_size]);
G = computeGeometry(G);

%% We need two rock structures for the dual-porosity model: one for the fracture
% system and one for the non-fractured rock matrix. Fractures usually will 
% have high permeability and low pore volume, while matrix usually will 
% have a higher porosity and decreased permeability. Therefore,
% interconnected fractures generally serves as a "highway" for fluid flow, 
% while the fluid transfer with the matrix will happen at a larger
% timescale
kf = 10000*darcy;
phif = 0.2;

km = 1*darcy;
phim = 0.5;

% Rock Structures
rock_fracture = makeRock(G, kf, phif);
rock_matrix = makeRock(G, km, phim);

%% Pressures of injector and producer wells
pprod = 0;
pinj = 1*barsa;

%% We also need two fluid structures for the dual-porosity model. Fractures are
% usually set as having zero capillary pressure, but we don't define any 
fluid_matrix = initSimpleADIFluid( 'mu' , [1, 2, 0] .* centi*poise, ...
                           'rho', [1000, 700, 0] .* kilogram/meter^3, ...
                           'n'  , [2, 2, 0]);
fluid_fracture = fluid_matrix;

bp = 1;
bi = 1.01;

m = (bi-bp)/(pinj-pprod);
b = bp-m*pprod;

fluid_fracture.bW = @(p)(b+m*p);
fluid_fracture.bO = @(p)(b+m*p);


%% Define a first-order transfer-function between fracture and matrix.
% First-order transfer functions are popular because they proportional to 
% a potential difference(in this case, the saturation difference, but it 
% could be pressure difference, capillary pressure difference, etc.).
% The function handles are stored in the dp_info struct with the names 
% transfer_water and transfer_oil, and receive (pwf,pwm,swf,swm,pof,pom,sof,som) 
% as input arguments. Pressures and saturations in the fracture and matrix 
% will drive the fluid exchange between the two continua. Here we define
% the transfer of oil phase as the opposite of the transfer of water phase,
% meaning that any amount of water that we transfer from fracture to the
% matrix will displace a mass of oil that will be transferred to the fractures 
beta = 5e-05;
dp_info.transfer_water = @(pwf,pwm,swf,swm,...
                       pof,pom,sof,som) ( beta * (swf-swm));

dp_info.transfer_oil = @(pwf,pwm,swf,swm,...
                     pof,pom,sof,som) ( -beta * (swf-swm));


%% Add boundary conditions
pv = poreVolume(G,rock_matrix);
bc = [];
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

%bc = addBC(bc, left, 'pressure', pinj, 'sat', [1 0]);
bc = addBC(bc, left, 'flux', sum(pv)/(10*day), 'sat', [1 0])
bc = addBC(bc, right, 'pressure', pprod, 'sat', [0 1]);


%% Create the model. TwoPhaseOilWaterDPModel is a modified version of 
% TwoPhaseOilWaterModel that adds dual porosity behaviour
gravity off

model = TwoPhaseOilWaterDPModel(G, rock_fracture, fluid_fracture,...
                            rock_matrix, fluid_matrix, dp_info);

%% Initialize the field with a constant pressure and fully saturated by oil
state0.pressure = ones(G.cells.num,1)*pprod;
state0.s = repmat([0 1],G.cells.num,1);
state0.swm = zeros(G.cells.num,1);
state0.pom = ones(G.cells.num,1)*pprod;

%% Initialize the well solution
state0.wellSol= initWellSolAD([], model, state0);
state = state0;

solver = NonLinearSolver();

%% Handles to pictures that are going to be plotted at each timestep
fig1 = figure('Position',[100,0,600,800]);
fig2 = figure('Position',[700,0,600,800]);

%% Simulate the models
Time = 10* day;
n = 5;
dT = Time/n;
for i = 1:n
    %% Advancing fields
    state = solver.solveTimestep(state, dT, model, 'bc', bc);
    disp(i)

    %% Plotting fields
    figure(fig1)
    subplot(2,1,1)
    title(['Saturation of Oil in the Fractures @ Time = ' num2str(i*dT/day) ' days'])
    p = plotCellData(G,1-state.s(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar
    caxis([0 1])

    figure(fig1)
    subplot(2,1,2)
    title(['Saturation of Oil in the Matrix @ Time = ' num2str(i*dT/day) ' days'])
    p = plotCellData(G,1-state.swm(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar
    caxis([0 1])

    figure(fig2)
    subplot(2,1,1)
    title(['Transfer Rate of Water @ Time = ' num2str(i*dT/day) ' days'])
    p = plotCellData(G,state.Twm(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar

    figure(fig2)
    subplot(2,1,2)
    title(['Transfer Rate of Oil @ Time = ' num2str(i*dT/day) ' days'])
    p = plotCellData(G,state.Tom(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar

    drawnow;
end
