function [] = DPvsMsRSB_2ph()  

    clc;
    close all;

    %% Read the problem from a deckfile
    mrstModule add ad-props ad-core ad-blackoil ad-fi

    %% Create grid 
    cells = 39;
    celldim = [4*cells cells];
    physdim = [4*cells cells];
    G = cartGrid(celldim, physdim);
    G = computeGeometry(G);
    
    %% Create rock
    % Base porosity and permeability
    kf = 10000*darcy;
    phif = 0.5;
    km = 1*darcy;
    phim = 0.5;
    
    % Two layer grid
    rock_fracture = makeRock(G, kf, phif);
    rock_matrix = makeRock(G, km, phim);
    
    %% Define pressures (drawdown)
    pres = 100*barsa;
    pinj = 120*barsa;

    %% Create fluid
    % Same fluid properties for both continua
    fluid_fracture = initSimpleADIFluid( 'mu' , [1, 2, 0] .* centi*poise, ...
        'rho', [1000, 700, 0] .* kilogram/meter^3, ...
        'n'  , [2, 2, 0]);
    Bwb = 1;
    Bwi = 1.001;
    
    Bob = 1;
    Boi = 1.001;
    
    fluid_fracture.bW = @(p) equation_line( p,[pres,Bwb],[pinj,Bwi] );
    fluid_fracture.bO = @(p) equation_line( p,[pres,Bob],[pinj,Boi] );
    
    fluid_matrix = fluid_fracture;
    
        %% Add wells
    % We will include two wells, both the injector and the producer controlled
    % by bottom-hole pressure. The injector and producer are located at the
    % mid-point of the left and right boundary of the domain.

    W = addWell([], G, rock_matrix, (floor(cells/2)-1)*celldim(1)+1 ,'InnerProduct', 'ip_tpf','Type', ...
        'bhp', 'Val', 120*barsa,'Radius', .1, 'Comp_i', [1, 0]);
    W = addWell(W, G, rock_matrix, floor(cells/2)*celldim(1), 'InnerProduct', 'ip_tpf', 'Type', ...
        'bhp' , 'Val', pres, 'Radius', .1, 'Comp_i', [0, 1]);
% 
%     %% Set Boundary Conditions
%     bc_faces = find(abs(G.faces.centroids(:,1)-0)<=0.000001 &...
%                     G.faces.centroids(:,3)>=z_size-dz & ...
%                     G.faces.centroids(:,3)<=z_size);
%     
%     bc = [];
%     bc = addBC(bc,bc_faces,'pressure',pres,'sat',[0,1]);
%     bc = pside(bc,G,'ZMax',0,'sat',[1,0]);

    %% Schedule
    Time = 10*day;
    nt = 75;
    dt = repmat(1/nt*Time, nt, 1);
    schedule = simpleSchedule(dt, 'W', W);
    
    % Set up reservoir
    gravity on
    
    %% Create the models
    dp_info.shape_factor = 1e-09;
    factor = 1e-06;
    dp_info.transfer_water = @(pf,pm,sf,sm)transfer_water(pf,pm,sf,sm,factor);
    dp_info.transfer_oil = @(pf,pm,sf,sm)transfer_oil(pf,pm,sf,sm,factor);
    model = TwoPhaseOilWaterDPModel(G, rock_fracture, fluid_fracture, rock_matrix, fluid_matrix, dp_info);
    
    %% Initializing state
     state.pressure = ones(G.cells.num,1)*0;
     state.s = repmat([0 1],G.cells.num,1);
     
     state.swm = state.s;
     state.pom = state.pressure;
    
    %% Simulate the models
    [wellSols, states] = simulateScheduleAD(state, model, schedule);
    
    %% Plotting Information
    ti = 0;
    t = zeros(length(states),1);
    ps = zeros(length(states),1);
%     fig1 = figure('Position',[100,100,600,600]);
%     fig2 = figure('Position',[700,100,600,600]);
%     fig3 = figure('Position',[1300,100,600,600]);
    for i = 1:1:length(states)
        ti = ti + dt(i);
        t(i) = ti;
        ps(i,1)= states{i,1}.s(W(2).cells,1);
    end
    
    plotNo = 1;
    figure; hold on; colormap(flipud(winter))
 
    for j = nt/3:nt/3:nt
        % Plot saturation
        r = 0.01;
        subplot(1,3,plotNo);
        plotCellData(G,states{j,1}.s(:,1),'EdgeColor','none');
        title(['S_w after: ', num2str(t(j)/day), 'days']);
        view(0,90); caxis([0 1]);
        axis tight;
        plotNo = plotNo+1;
        text(G.cells.centroids(W(1).cells,1),G.cells.centroids(W(1).cells,2),'x')
        text(G.cells.centroids(W(2).cells,1),G.cells.centroids(W(2).cells,2),'o')
    end
        

    
    figure;
    title(['Saturation at producer']);
    plot(t(:,1)/day, ps(:,1),'-x');
    ylabel('Saturation at producer');
    xlabel('t [days]'); 
    %set(gca,'XGrid','on','YGrid','on');
    axis([0 Time/day 0 1]);

end

function [T] = transfer_water(pwf,pwm,swf,swm,factor)

    T = (pwf-pwm)*0;
    
end

function [T] = transfer_oil(pof,pom,sof,som,factor)
    
    satmin = 0.1;
    T = factor * ...
        double(sof > satmin) .* ...
        (sof-som);
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


