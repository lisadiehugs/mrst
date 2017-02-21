%% Hadi with boundary conditions
% we apply pressure drop boundary conditions to a quasi 1D domain. This 
% test is for comparison with the hfm/ F-MsRSB approach
clc;
close all; clear all;

%% Include the essential modules. dual-porosity model has the functionality
% we need.
mrstModule add ad-props ad-core ad-blackoil ad-fi blackoil-sequential dual-porosity

%% Create the grid: a simple cartesian grid, 100m x 100m and 50 blocks on x
% and y directions.

% Two layer grid
y_size = 50;
x_size = 100;

Ny = 1; % 1D problem for dual-porosity
Nx = 20;

dx = x_size/Nx;


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
phif = 0.5;

km = 1e-15*darcy;
phim = 0.5;

% Rock Structures
rock_fracture = makeRock(G, kf, phif);
rock_matrix = makeRock(G, km, phim);

%% Pressures of injector and producer wells
pprod = 1*barsa;
pinj = 10*barsa;

%% We also need two fluid structures for the dual-porosity model. Fractures are
% usually set as having zero capillary pressure, but we don't define any 
fluid_matrix = initSimpleADIFluid( 'mu' , [1, 2, 0] .* centi*poise, ...
                           'rho', [1000, 700, 0] .* kilogram/meter^3, ...
                           'n'  , [1, 1, 0]);
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
beta = 5e-06;
dp_info.transfer_water = @(pwf,pwm,swf,swm,...
                       pof,pom,sof,som)transfer_water(pwf,pwm,swf,swm,...
                       pof,pom,sof,som,beta,dx);
dp_info.transfer_oil = @(pwf,pwm,swf,swm,...
                       pof,pom,sof,som)transfer_oil(pwf,pwm,swf,swm,...
                       pof,pom,sof,som,beta,dx);
% dp_info.transfer_water = @(pwf,pwm,swf,swm,...
%                        pof,pom,sof,som) ( beta * (swf-swm));
% 
% dp_info.transfer_oil = @(pwf,pwm,swf,swm,...
%                      pof,pom,sof,som) ( -beta * (swf-swm));


%% Add boundary conditions
pv = poreVolume(G,rock_matrix);
bc = [];
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = addBC(bc, left, 'pressure', pinj, 'sat', [1 0]);
%bc = addBC(bc, left, 'flux', sum(pv)/(10*day), 'sat', [1 0])
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
fig3 = figure('Position',[700,0,600,800]);

%% Simulate the models
Time = 0.03* day;
n =250;
dT = Time/n;
t = 0;
% In-/Outflow face
xf = G.faces.centroids(:,1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);
totfluxn = zeros(n,1);
totfluxw = zeros(n,1);
sout = zeros(n,1);
for i = 1:n
    %% Advancing fields
    state = solver.solveTimestep(state, dT, model, 'bc', bc);
    disp(i)
    
    % get outflow
    totfluxn(i) = state.flux(right,2);
    totfluxw(i) = state.flux(right,1);
    sout(i) = state.s(end,1);

    %% Plotting fields
%     if (mod(i,5)~=0)
%         figure(fig1)
%         subplot(2,1,1)
%         title(['Saturation of Oil in the Fractures @ Time = ' num2str(i*dT/day) ' days'])
%         p = plotCellData(G,1-state.s(:,1));
%         p.EdgeAlpha = 0;
%         axis equal tight off
%         view(-20, 50)
%         colorbar
%         caxis([0 1])
%         
%         figure(fig1)
%         subplot(2,1,2)
%         title(['Saturation of Oil in the Matrix @ Time = ' num2str(i*dT/day) ' days'])
%         p = plotCellData(G,1-state.swm(:,1));
%         p.EdgeAlpha = 0;
%         axis equal tight off
%         view(-20, 50)
%         colorbar
%         caxis([0 1])
%         
%         drawnow;
%     end
    
%     figure(fig2)
%     plot(t/day,state.s(end,1),'ro')
%     hold on
%     grid on
    
    figure(fig3)
    title(['Pressure of Oil in the Matrix @ Time = ' num2str(i*dT/day) ' days'])
    p = plotCellData(G,state.pressure);
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar

    drawnow;
    
%     pause(3)
    
    t=t+dT;
    
end

 %% plot saturation and flux at outflow
    figure(7);
    subplot(1,2,1)
    plot(linspace(0, Time, n)/day, sout,'k-o');
    leg = legend('Fine-scale','Location','Best');
    ylabel('Saturation at Outflow boundary [-]');
    xlabel('Time [days]');
    %saveas(7, ['watercutKf=' num2str(K_frac) 'Km=' num2str(K_m) '.epsc']);
    subplot(1,2,2);
    plot(linspace(0, Time, n)/day, totfluxw.*(-1),'-x', linspace(0, Time, n)/day, totfluxn.*(-1),'-o');
    leg = legend('water', 'oil','Location','Best');
    ylabel('Flux at Outflow boundary [-]');
    xlabel('Time [days]');

function Tw = transfer_water(pwf,pwm,swf,swm,pof,pom,sof,som,beta,dx) 
    if (isprop(pwf,'val'))
        dpwf = diff(pwf.val)/dx;
    else
        dpwf = diff(pwf)/dx;
    end
    dpwf = [dpwf(1);dpwf];
    %disp(dpwf)
    Tw = beta * (swf-swm);
    
end

function To = transfer_oil(pwf,pwm,swf,swm,pof,pom,sof,som,beta,dx) 

    Tw = transfer_water(pwf,pwm,swf,swm,pof,pom,sof,som,beta,dx) ;
    To = -Tw;
    
end
