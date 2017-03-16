%% F-MsRSB applied to a heterogeneous 2D domain
% Two-phase example modeling water injection through a quarter 5-spot into a
% 2-dimensional fractured porous media. The flow problem is solved both by a
% fine-scale and a multiscale solver.

% Notice that you need to have Metis installed to get this example to work.
% To get Metis working, you also need to set the global variable METISPATH.
% This can be done in your 'startup_user.m' file.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Call run2phFiveSpot
K_frac = [1e3;5e3;1e4;5e4; 1e5;1e6];
basis_ms = cell(size(K_frac, 1)+1, 1);
for i = 1:size(K_frac, 1)
    [G, CG, basis, ~] = runDPVglBasis_fitDP_fluxBC(K_frac, i);
    basis_ms{i+1} = basis{1};
    if (i == 1)
        basis_ms{1}=basis{2};
    end
end
K_frac = [0;K_frac];
plotCompareChooseBasisTimestep(basis_ms, G, CG, K_frac, 'K_frac');



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