function [G, CG, basis, pv_out] = run2phBC(K_frac, curK)
    output = true;
    if (size(K_frac, 1) > 2)
        output = false;
    end
    %% Grid and fracture line(s)
    % Construct a Cartesian grid comprising 100-by-100 cells, where each cell
    % has dimension 1-by-1 m^2. Define 1 fracture line, roughly 85 m in length,
    % extending diagonally through the centre of the domain.

    celldim = [90 90];
    physdim = [9 9];
    G = cartGrid(celldim, physdim);
    G = computeGeometry(G);

    fl = [  % // to x-axis
         0, 3, 9, 3;
         0, 6, 9, 6;
         % // to y-axis
         4.5 , 0 , 4.5, 9];

    %% Process fracture line(s)
    % Using the input fracture lines, we identify independent fracture networks
    % comprising of connected lines. In this example, there is only 1 fracture
    % network consisting of 1 fracture line. We also identify the fine-cells
    % in the matrix containing these fractures. Fracture aperture is set to
    % 0.04 meters. The matrix grid and fracture lines are plotted.

    dispif(mrstVerbose, 'Processing user input...\n\n');
    [G,fracture] = processFracture2D(G,fl);
    fracture.aperture = 1/25; % Fracture aperture

    %% Compute CI and construct fracture grid
    % For each matrix block containing a fracture, we compute a fracture-matrix
    % conductivity index (CI) in accordance with the hierarchical fracture
    % model (a.k.a. embedded discrete fracture model). Following that, we
    % compute a fracture grid where the fracture cell size is defined to be
    % 1 m. Next, the fracture grid is plotted on top of the matrix grid.

    dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
    G = CIcalculator2D(G,fracture);
    min_size = 0.4; cell_size = 0.5; % minimum and average cell size.
    [G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);

    %% Set rock properties in fracture and matrix
    % Set the permeability (K) as 1 Darcy in the matrix and 10000 Darcy in the
    % fractures. Additionally, set the porosity of the matrix and fractures to
    % 50%.

    dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
    G.rock.perm = ones(G.cells.num,1)*darcy;
    G.rock.poro = 0.5*ones(G.cells.num, 1);
    poro_frac = 0.2;
    G = makeRockFrac(G, K_frac(curK), 'permtype','homogeneous','porosity',poro_frac);

    %% Define fluid properties
    % Define a two-phase fluid model without capillarity. The fluid model has
    % the values:
    %
    % * densities: [rho_w, rho_o] = [1000 700] kg/m^3
    % * viscosities: [mu_w, mu_o] = [1 2] cP.
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

    %% Add BC
    % Set boundary condition P = 10 bars on the left face and P = 1 bar on the
    % right face of the domain.

    bc = [];
    xf = G.Matrix.faces.centroids(:, 1);
    left = find(abs(xf - min(xf)) < 1e-4);
    right = find(abs(xf - max(xf)) < 1e-4);

    bc = addBC(bc, left, 'pressure', 10*barsa, 'sat', [1 0]);
    bc = addBC(bc, right, 'pressure', 1*barsa, 'sat', [0 1]);

    boundfaces=findfracboundaryfaces2d(G,1e-5);
    bc = addBC(bc,boundfaces.West,'pressure',10*barsa, 'sat', [1 0]);
    bc = addBC(bc,boundfaces.East,'pressure',1*barsa, 'sat', [0 1]);

    %% Initialize state variables
    % Once the transmissibilities are computed, we can generate the
    % transmissiblity matrix 'A' using the 'two-point flux approximation'
    % scheme and initialize the solution structure.

    dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
    state  = initResSol (G, 1*barsa, [0 1]);
    [A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true); 

    %% Setup multiscale grids 
    % Next, we define a 10-by-10 matrix coarse grid and 2 coarse blocks (or
    % coarse degrees of freedom) in the fracture. Additionally, we also define
    % the support regions for the fracture and matrix basis functions. The
    % matrix and fracture coarse grids are plotted.

    dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');
    coarseDims = [9 9]; % coarsening factor in each direction
    dof_frac = 18; % Fracture dof per fracture network
    [CG, CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims,...
        'dof_frac',dof_frac, ...
        'coarseNodeOption', 'useCoarseCellEndPoints',...
        'fracSupportRadius', 6); % 'paddedPartition', true, 
    figure(1); plotFractureCoarseGrid2D(G,CG.partition,F)

    %% Compute initial pressure

    dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
    state_fs = incompTPFA(state, G, T, fluid,  ...
        'bc',bc, 'MatrixOutput', true, 'use_trans',true);

    %% Compute basis functions
    % Using the transmissibility matrix 'A' we compute the basis functions for
    % each fracture and matrix coarse block using the restriction smoothed
    % basis method. Note that the matrix 'A' does not contain any source terms
    % or boundary conditions. They are added to the coarse linear system when
    % computing the multiscale pressure in the next section.

    dispif(mrstVerbose, 'Computing basis functions...\n\n');
    basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
    if (output)
        figure; plotToolbar(G,basis_sb.B,'filterzero',true);
        prm = log10(G.rock.perm); mx = max(prm); mn = min(prm);
        G.nodes.z = 1e-3*ones(G.nodes.num,1);
        plotCellData(G,(prm-mn)./(mx-mn),'EdgeColor','none','FaceAlpha',.4);
        G.nodes = rmfield(G.nodes,'z');
        plotGrid(CG,'FaceColor','none');
        axis tight; colorbar;
        title('Basis functions plotted in the matrix');
    end

    %% Compute multiscale solution

    dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
    [state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
        'bc', bc,'use_trans',true);

    %% Incompressible Two-Phase Flow
    % We solve the two-phase system using a sequential splitting in which the
    % pressure and fluxes are computed by solving the flow equation and then
    % held fixed as the saturation is advanced according to the transport
    % equation. Note that the flow equations are solved in the fine scale as
    % well as on the coarse scale using an algebraic multiscale strategy. The
    % multiscale flux field obtained at fine scale resolution is reconstructed
    % to be conservative before solving the transport equation. This procedure
    % is repeated for a given number of time steps (here we use 45 equally
    % spaced time steps). The error introduced by this splitting of flow and
    % transport can be reduced by iterating each time step until e.g., the
    % residual is below a certain user-prescribed threshold (this is not done
    % herein).

    pv     = poreVolume(G,G.rock);
    nt     = 50;
    t200    = 0.7*(sum(pv)/sum(state_fs.flux(left)));
    Time   = t200;
    dT     = Time/nt;
    dTplot = Time/3;
    N      = fix(Time/dTplot);

    pvi = zeros(nt,1);
    sol_fs = cell(nt,1); sol_ms = cell(nt,1);
    e = zeros(nt,1); pms = zeros(nt,1); pfs = zeros(nt,1);
    
    plotNo = 1; hfs = 'Reference: '; hms = 'F-MsRSB: ';
    t  = 0; 
    B = basis_sb.B;
    R = controlVolumeRestriction(CG.partition);
    count = 1;
    hwb = waitbar(0,'Time loop');
    while t < Time
        state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'bc', bc, 'Trans', T,'verbose',true);
        state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'bc', bc, 'Trans', T);
        % Check for inconsistent saturations
        s = [state_fs.s(:,1); state_ms.s(:,1)];
        assert(max(s) < 1+eps && min(s) > -eps);

        % Update solution of pressure equation.
        state_fs  = incompTPFA(state_fs , G, T, fluid, 'bc', bc, 'use_trans',true);

        %-------------------------------Multiscale----------------------------%

        A = getSystemIncompTPFA(state_ms, G, T, fluid, 'use_trans', true);
        B = iteratedJacobiBasis(A, CG, 'interpolator', B); 
        basis_sb = struct('B', B, 'R', R);
        state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb, 'bc', bc,'use_trans',true);

        %---------------------------------------------------------------------%

        sol_fs{count,1} = state_fs; sol_ms{count,1} = state_ms;
        % Increase time
        t = t + dT;
        waitbar(t/Time,hwb);
        count = count + 1;
        if ( t < plotNo*dTplot && t < Time), continue, end
    
    % Plot saturation
    heading = [num2str(convertTo(t,day)),  ' days'];
    
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.44]), cla
    plotCellData(G, state_fs.s(:,1), 'edgealpha', 0.1);
    colormap(flipud(jet))
    view(0,90), axis equal off, title([hfs heading],'FontSize',8)
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.44]), cla
    plotCellData(G, state_ms.s(:,1), 'edgealpha', 0.1);
    colormap(flipud(jet))
    view(0,90), axis equal off, title([hms heading],'FontSize',8)
    
    plotNo = plotNo+1;
    drawnow

    end
    close(hwb);

    %% Plot saturations
    if (output)
        OutflowfacesMat = right;
        OutflowfacesFrac = boundfaces.East';
        Outflowfaces = [OutflowfacesMat;OutflowfacesFrac];

        figure;
        set(gcf,'Position',[200 200 1000 450]);
        subplot(1,2,1)
        plot(G.faces.centroids(OutflowfacesMat,2), state_fs.flux(OutflowfacesMat),'-o', ...
            G.faces.centroids(OutflowfacesMat,2), state_ms.flux(OutflowfacesMat),'-*');
        leg = legend('Fine-scale','Multiscale','Location','Best');
        ylabel('Flux at Outflow boundary [m^3/s]');
        xlabel('y [m]'); 
        title(['Matrix cells'])
        set(gca,'XGrid','on','YGrid','on');
        axis([0 physdim(2) 0 1.1*max([state_fs.flux(OutflowfacesMat); state_ms.flux(OutflowfacesMat)])])
        subplot(1,2,2)
        plot(G.faces.centroids(OutflowfacesFrac,2), state_fs.flux(OutflowfacesFrac),'o',G.faces.centroids(OutflowfacesFrac,2), state_ms.flux(OutflowfacesFrac),'*');
        leg = legend('Fine-scale','Multiscale','Location','Best');
        ylabel('Flux at Outflow boundary [m^3/s]');
        xlabel('y [m]'); 
        title(['Fracture Cells'])
        set(gca,'XGrid','on','YGrid','on');
        axis([0 physdim(2) 0 1.1*max([state_fs.flux(OutflowfacesFrac); state_ms.flux(OutflowfacesFrac)])])
    end

%     %% Plot water saturation at producer 
% 
%     figure;
%     plot(pvi,pfs(:,1),'-o',pvi,pms(:,1),'--*');
%     leg = legend('Fine-scale','Multiscale','Location','Best');
%     ylabel('Saturation at producer');
%     xlabel('PVI [%]'); 
%     set(gca,'XGrid','on','YGrid','on');
%     axis tight

%     %% Plot error in saturation 
%     if (output)
%         figure;
%         plot(pvi,e*100, '--+b');
%         ylabel('e [%]')
%         xlabel('PVI [%]'); 
%         set(gca,'XGrid','on','YGrid','on');
%         axis tight
% 
%         e_eq = '$$ e = \frac{ \sum ( |S_w^{fs}-S_w^{f-msrsb}| \times pv) }{ \sum (S_w^{fs} \times pv) } $$';
%         title(e_eq,'interpreter','latex');
%     end
        basis = basis_sb;
        pv_out = [pvi,pfs,pms];
end
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
