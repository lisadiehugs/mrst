function [G, CG, basis_out, pvOut] = runDPVglBasis_fitDP_fluxBC_II(K_frac, curK)
%% test case for comparison with dual-porosity model

    %% Grid and fracture line(s)
    % Construct a Cartesian grid comprising 200-by-100 cells, where each cell
    % has dimension 1-by-1 m^2. Define crossing pattern of fracture lines.

    celldim = [200 100];
    physdim = [100 50];
    G = cartGrid(celldim, physdim);
    G = computeGeometry(G);
    dx = 100/120;

    fl = [0, 25, 25, 50;
          75, 0, 100, 25;
          %l.o.r.u.
          0,50,25,0;
          25,50,50,0;
          50,50,75,0;
          75,50,100,0;

          % "inner" fractures 
          0,0,50,50;
          25, 0, 75, 50;
          50, 0, 100, 50];

    %% Process fracture line(s)
    % Using the input fracture lines, we identify independent fracture networks
    % comprising of connected lines. In this example, there is only 1 fracture
    % network consisting of 1 fracture line. We also identify the fine-cells
    % in the matrix containing these fractures. Fracture aperture is set to
    % 0.04 meters. The matrix grid and fracture lines are plotted.

    dispif(mrstVerbose, 'Processing user input...\n\n');
    [G,fracture] = processFracture2D(G,fl);
    fracture.aperture = 1/25; % Fracture aperture
    gravity off

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
    
    K_m = 1;
    dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
    G.rock.perm = K_m*ones(G.cells.num,1)*darcy;
    G.rock.poro = 0.2*ones(G.cells.num, 1); 
    poro_frac = 0.5;
    G = makeRockFrac(G, K_frac(curK), 'permtype','homogeneous','porosity', poro_frac);

    %% Define fluid properties
    % Define a two-phase fluid model without capillarity. The fluid model has
    % the values:
    %
    % * densities: [rho_w, rho_o] = [1000 700] kg/m^3
    % * viscosities: [mu_w, mu_o] = [1 2] cP.
    % * corey-coefficient: [2, 2] = [2 2].

    fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
        'rho', [1000, 700] .* kilogram/meter^3, ...
        'n'  , [  1,   1]);

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
    boundfaces=findfracboundaryfaces2d(G,1e-5);
    Inflowfaces = [left;boundfaces.West'];
    flux = 1e3/(1*day);
    InflowFrac = flux/size(boundfaces.East',1); %G.faces.areas(boundfaces.West(1))*flux/sum(G.faces.areas(Inflowfaces),1);
    InflowMat = G.faces.areas(left(1))*flux/sum(G.faces.areas(Inflowfaces),1);


    %bc = addBC(bc, left, 'flux', InflowMat, 'sat', [1 0]);
    bc = addBC(bc, right, 'pressure', 1*barsa, 'sat', [0 1]);

    bc = addBC(bc, boundfaces.West, 'flux', InflowFrac, 'sat', [1 0]);
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
    coarseDims = [9 5]; % coarsening factor in each direction
    p = partitionMatrix(G, 'coarseDims' , coarseDims);
    G_boundcells = nonzeros(G.Matrix.faces.neighbors(find(~all(G.Matrix.faces.neighbors,2)),:));
    CG_boundcells = unique(p(G_boundcells));
    CC = [1:coarseDims(1)*coarseDims(2)]';
    CG_innercells = find(~ismember(CC, CG_boundcells));
    leftCC = CC(0< mod(1:coarseDims(1)*coarseDims(2),coarseDims(1)) & mod(1:coarseDims(1)*coarseDims(2),coarseDims(1)) <5);
    rightCC = CC( mod(1:coarseDims(1)*coarseDims(2),coarseDims(1)) < 1 | mod(1:coarseDims(1)*coarseDims(2),coarseDims(1)) >5);
    p(find(ismember(p, CG_innercells) & ismember(p, leftCC))) = coarseDims(1)*coarseDims(2)+1;
    p(find(ismember(p, CG_innercells) & ismember(p, rightCC))) = coarseDims(1)*coarseDims(2)+2;
    [~,~,pm] = unique(p);
    dof_frac = 32; % Fracture dof per fracture network
    [CG, CGf] = getRsbGridsHFM(G, fracture.network, 'pm' , pm,...
        'dof_frac',dof_frac, ...
        'fracSupportRadius', 8, 'removeCenters', true); % 
    if (curK ==1)
        state_MsRSB  = initResSol (G.Matrix, 1*barsa, [0 1]);
        p = partitionUI(G.Matrix, coarseDims);
        T_MsRSB = computeTrans(G.Matrix, G.Matrix.rock);
        CG_MsRSB = coarsenGeometry(generateCoarseGrid(G.Matrix, p));
        A_MsRSB = getSystemIncompTPFA(state_MsRSB, G.Matrix, T_MsRSB, fluid);
        basis_MsRSBIni = getMultiscaleBasis(CG_MsRSB, A_MsRSB, 'type', 'rsb', 'useMex', false, 'Iterations', 1000);
    end
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
%         figure; plotToolbar(G,basis_sb.B,'filterzero',true);
%         prm = log10(G.rock.perm); mx = max(prm); mn = min(prm);
%         G.nodes.z = 1e-3*ones(G.nodes.num,1);
%         plotCellData(G,(prm-mn)./(mx-mn),'EdgeColor','none','FaceAlpha',.4);
%         G.nodes = rmfield(G.nodes,'z');
%         plotGrid(CG,'FaceColor','none');
%         axis tight; colorbar;
%         title('Basis functions plotted in the matrix');

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

    nt     = 1;
    t200    = 0.01*day; 
    Time   = t200;
    dT     = Time/nt;
    dTplot = Time/3;
    N      = fix(Time/dTplot);

    pvi = zeros(nt,1);
    s_fsout = zeros(nt,1); % stores accumulated saturation along outflow boundary
    s_msout = zeros(nt,1);
    s_fs = cell(nt,1); % stores accumulated saturation for domain
    s_ms = cell(nt,1);
    sfrac_fs = cell(nt,1);
    sfrac_ms = cell(nt,1);

    sol_fs = cell(nt,1); sol_ms = cell(nt,1);
    basis = cell(nt,1);
    basis{1,1} = basis_sb;
    if (curK==1)
        basis_MsRSB{1,1} = basis_MsRSBIni;
    end
    ef = zeros(nt,1); pms = zeros(nt,1); pfs = zeros(nt,1);
    plotNo = 1; hfs = 'Reference: '; hms = 'F-MsRSB: ';
    t  = 0; 
    B = basis{1,1}.B;
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
        basis{count,1} = struct('B', B, 'R', R);
        state_ms = incompMultiscale(state_ms, CG, T, fluid, basis{count,1}, 'bc', bc,'use_trans',true);
        if (curK==1)
            basis_MsRSB{count,1} = getMultiscaleBasis(CG_MsRSB, A_MsRSB, 'type', 'rsb', 'useMex', false, 'Iterations', 1000);
        end
        %---------------------------------------------------------------------%

        sol_fs{count,1} = state_fs; sol_ms{count,1} = state_ms;
        

        % saturation in domain for matrix (accumulate values into nmx blocks)
        nmx = 8;
        em = 0:physdim(1)/nmx:physdim(1);
        ym = discretize(G.cells.centroids(1:G.Matrix.cells.num, 1),em);
        s_fs{count} = accumarray(ym, state_fs.s(1:G.Matrix.cells.num,1))./histcounts(G.cells.centroids(1:G.Matrix.cells.num, 1),em)';
        s_ms{count} = accumarray(ym, state_ms.s(1:G.Matrix.cells.num,1))./histcounts(G.cells.centroids(1:G.Matrix.cells.num, 1),em)';
        
        % saturation in domain for fracture (accumulate values into nfx blocks)
        nfx = 20;
        ef = 0:physdim(1)/nfx:physdim(1);
        yf=discretize(G.cells.centroids((G.Matrix.cells.num+1):end, 1), ef);
        sfrac_fs{count} = accumarray(yf, state_fs.s(G.Matrix.cells.num+1:end,1))./ histcounts(G.cells.centroids(G.Matrix.cells.num+1:end, 1),ef)';
        sfrac_ms{count} = accumarray(yf, state_fs.s(G.Matrix.cells.num+1:end,1))./ histcounts(G.cells.centroids(G.Matrix.cells.num+1:end, 1),ef)';
        % integrate and store saturation at outflow boundary
        s_fsout(count) = s_fs{count}(end);
        s_msout(count) = s_ms{count}(end);
%         if (mod(count,2)~=0)
%             figure(4)
%             set(4, 'Position', [100, 100, 900, 400])
%             subplot(1,2,1)
%             plot(em(1:end-1)+em(2)/2, s_fs{count}, '-o', em(1:end-1)+em(2)/2, s_ms{count}, '-x');
%             ylabel('Saturation [-]');
%             xlabel('x [m]');
%             legend('Fine-scale','Multiscale','Location','Best');
%             axis([0 physdim(1) 0 1])
%             title(['Saturation in Matrix (' num2str(nmx) ' cells) - t = ' num2str(t/day) ' days'])
%             subplot(1,2,2)
%             plot(ef(1:end-1)+ef(2)/2, sfrac_fs{count}, '-o', ef(1:end-1)+ef(2)/2, sfrac_ms{count}, '-x');
%             legend('Fine-scale','Multiscale','Location','Best');
%             ylabel('Saturation [-]');
%             xlabel('x [m]');
%             axis([0 physdim(1) 0 1])
%             title(['Saturation in Fractures (' num2str(nfx) ' cells) - t = ' num2str(t/day) ' days'])
%         end
        % Increase time
        t = t + dT;
        waitbar(t/Time,hwb);
        count = count + 1;
        if ( t < plotNo*dTplot && t < Time), continue, end
%         saveas(4, ['satFrontKf=' num2str(K_frac) 'Km=' num2str(K_m) '_' num2str(plotNo) '.emf']);
%         figure(5);
%         set(5, 'Position', [100, 100, 1100, 500])
%         % Plot saturation
%         heading = [num2str(convertTo(t,day)),  ' days'];
% 
%         r = 0.01;
%         subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.44]), cla
%         plotCellData(G, state_fs.s(:,1), 'edgealpha', 0.1);
%         colormap(flipud(jet))
%         view(0,90), axis equal off, title([hfs heading],'FontSize',8)
% 
%         subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.44]), cla
%         plotCellData(G, state_ms.s(:,1), 'edgealpha', 0.1);
%         colormap(flipud(jet))
%         view(0,90), axis equal off, title([hms heading],'FontSize',8)
% 
%         plotNo = plotNo+1;
%         drawnow

    end
%     saveas(5, ['satPlotKf=' num2str(K_frac) 'Km=' num2str(K_m) '.emf']);
    close(hwb);


%     %% Plot fluxes and error
%     OutflowfacesMat = right;
%     OutflowfacesFrac = boundfaces.East';
%     Outflowfaces = [OutflowfacesMat;OutflowfacesFrac];
% 
%     figure(6);
%     subplot(1,2,1)
%     plot(G.faces.centroids(OutflowfacesMat,2), state_fs.flux(OutflowfacesMat),'-o', ...
%         G.faces.centroids(OutflowfacesMat,2), state_ms.flux(OutflowfacesMat),'-*');
%     leg = legend('Fine-scale','Multiscale','Location','Best');
%     ylabel('Flux at Outflow boundary [m^3/s]');
%     xlabel('y [m]'); 
%     title(['Matrix cells'])
%     set(gca,'XGrid','on','YGrid','on');
%     axis([0 physdim(2) 0 1.1*max([state_fs.flux(OutflowfacesMat); state_ms.flux(OutflowfacesMat)])])
%     
%     subplot(1,2,2)
%     el2 = sqrt(sum(abs((state_fs.flux(Outflowfaces)-state_ms.flux(Outflowfaces)).*(state_fs.flux(Outflowfaces)-state_ms.flux(Outflowfaces))).*G.faces.areas(Outflowfaces),1)/(sum(G.faces.areas(Outflowfaces).*state_fs.flux(Outflowfaces).*state_fs.flux(Outflowfaces),1)));
%     set(gcf,'Position',[200 200 600 400]);
%     eflux = abs(state_fs.flux(OutflowfacesMat)-state_ms.flux(OutflowfacesMat))./state_fs.flux(OutflowfacesMat);
%     plot(G.faces.centroids(OutflowfacesMat,2),eflux, '--+b');
%     ylabel('e')
%     xlabel('y [m]'); 
%     set(gca,'XGrid','on','YGrid','on');
%     axis ([0 physdim(2) 0 1.1*max(eflux)])
% 
%     e_eq = ['$$ e = \frac{ \sum ( |f^{fs}-f^{f-msrsb}|) }{ \sum (f^{fs}) } \qquad e_f^{L^2} =$$ ' num2str(el2)];
%     title(e_eq,'interpreter','latex'); 
%     set(6, 'Position', [100, 100, 900, 400])
%     %saveas(6, ['errorKf=' num2str(K_frac) 'Km=' num2str(K_m) '.emf']);
%     %% plot saturation at outflow
%     figure(7);
%     plot(linspace(0, Time, nt)/day, s_fsout,'-o', linspace(0, Time, nt)/day, s_msout, '-*');
%     leg = legend('Fine-scale','Multiscale','Location','Best');
%     ylabel('Saturation at Outflow boundary [-]');
%     xlabel('Time [days]');
%    % saveas(7, ['watercutKf=' num2str(K_frac) 'Km=' num2str(K_m) '.emf']);
    
pvOut= [];
basis_out{1}=basis;
if (curK==1)
    basis_out{2} = basis_MsRSB;
end
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