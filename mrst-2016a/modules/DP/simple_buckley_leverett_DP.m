function [] = simple_buckley_leverett_DP()    
    
    clc;
    close all;

    %% Read the problem from a deckfile
    mrstModule add ad-props ad-core ad-blackoil ad-fi mrst-gui

    %% Create grid 

    % Two layer grid
    x_size = 100;
    y_size = 100;
    z_size = 30;

    Nx = 50;
    Ny = 1;
    Nz = 50;
    
    dx = x_size/Nx;
    dy = y_size/Ny;
    dz = z_size/Nz;

    G = cartGrid([Nx Ny Nz],[x_size y_size z_size]);
    G = computeGeometry(G);
    
    
    %% Create rock
    % Base porosity and permeability
    kf = 10000*milli*darcy;
    phif = 0.001;
    km = 0.0001*milli*darcy;
    phim = 0.4;
    
    % Two layer grid
    rock_fracture = makeRock(G, kf, phif);
    rock_matrix = makeRock(G, km, phim);
    
    %% Define pressures (drawdown)
    pres = 1000*psia;
    pinj = pres + 100*psia;
    
    W = addWell([], G, rock_matrix, 1,'InnerProduct', 'ip_tpf','Type', ...
        'bhp', 'Val', 120*barsa,'Radius', .1, 'Comp_i', [1, 0]);

    %% Create fluid
    % Same fluid properties for both continua
    fluid_fracture = initSimpleADIFluid('mu', [2 0.1], 'rho', [100 0.1 1], 'n', [2 2 1]);
    Bwb = 1;
    Bwi = 1.001;
    
    Bob = 1;
    Boi = 1.001;
    
    fluid_fracture.bW = @(p) equation_line( p,[pres,Bwb],[pinj,Bwi] );
    fluid_fracture.bO = @(p) equation_line( p,[pres,Bob],[pinj,Boi] );
    
    fluid_matrix = fluid_fracture;

    %% Set Boundary Conditions
    bc_faces = find(abs(G.faces.centroids(:,1)-0)<=0.000001 &...
                    G.faces.centroids(:,3)>=z_size-dz & ...
                    G.faces.centroids(:,3)<=z_size);
    
    bc = [];
    bc = addBC(bc,bc_faces,'pressure',pres,'sat',[0,1]);
    bc = pside(bc,G,'ZMax',0,'sat',[1,0]);

    %% Schedule
    dt = repmat(0.1*day, 100, 1);
    schedule = simpleSchedule(dt,'W', W, 'bc', bc);
    
    % Set up reservoir
    gravity on
    
    %% Create the models
    zero_transfer_layer = find(G.cells.centroids(:,3)<=2*z_size/3 &...
                                G.cells.centroids(:,3)>=z_size/3);
    dp_info.shape_factor = 1e-09;
    factor = 1e-06;
    dp_info.transfer_water = @(pf,pm,sf,sm)transfer_water(pf,pm,sf,sm,factor,zero_transfer_layer);
    dp_info.transfer_oil = @(pf,pm,sf,sm)transfer_oil(pf,pm,sf,sm,factor,zero_transfer_layer);
    model = TwoPhaseOilWaterDPModel(G, rock_fracture, fluid_fracture, rock_matrix, fluid_matrix, dp_info);
    
    %% Initializing state
    state.pressure = ones(G.cells.num,1)*pres;
    state.s = repmat([1 0],G.cells.num,1);
    
    state.swm = state.s;
    state.pom = state.pressure;
    
    %% Simulate the models
    [wellSols, states] = simulateScheduleAD(state, model, schedule);
    
    %% Plotting Information
    ti = 0;
    
    fig1 = figure('Position',[100,100,600,600]);
    fig2 = figure('Position',[700,100,600,600]);
    fig3 = figure('Position',[1300,100,600,600]);
    for i = 1:1:length(states)
        ti = ti + dt(i);
        
        figure(fig1)
        title(['Saturation @ T = ' num2str(ti/day) ' days'])
        p = plotCellData(G,1-states{i}.s(:,1));
        p.EdgeAlpha = 0;
        axis equal tight off
        view(-20, 10)
        colorbar
        caxis([0 1])
        
        figure(fig2)
        title(['Saturation @ T = ' num2str(ti/day) ' days'])
        p = plotCellData(G,1-states{i}.swm(:,1));
        p.EdgeAlpha = 0;
        axis equal tight off
        view(-20, 10)
        colorbar
        caxis([0 1])
        
        drawnow;
        
        figure(fig3)
        title(['Transfer @ T = ' num2str(ti/day) ' days'])
        p = plotCellData(G,states{i}.Tom(:,1));
        p.EdgeAlpha = 0;
        axis equal tight off
        view(-20, 10)
        colorbar
        
        drawnow;
        
        
    end

end

function [T] = transfer_water(pwf,pwm,swf,swm,factor,zero_transfer_layer)

    T = (pwf-pwm)*0;
    T(zero_transfer_layer) = 0;
    
end

function [T] = transfer_oil(pof,pom,sof,som,factor,zero_transfer_layer)
    
    satmin = 0.1;
    T = factor * ...
        double(sof > satmin) .* ...
        (sof-som);
    T(zero_transfer_layer) = 0;
end

function [ y ] = equation_line( x,p1,p2 )

    x1 = p1(1);
    y1 = p1(2);
    
    x2 = p2(1);
    y2 = p2(2);
    
    m = (y2-y1)/(x2-x1);
    
    b = y2-m*x2;
    
    y = m*x+b;

end


