function [info, present] = dataset_spe10()
% Info function for SPE10 dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    [info] = datasetInfoStruct(...
        'name', 'SPE10', ...
        'website', '', ...
        'fileurl', '', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', true, ...
        'cells',   1122000, ...
        'examples', {'spe10:SPE10_rock', ...
                     'spe10:SPE10_setup'}, ...
        'description', [ ...
            'The tenth SPE benchmark was initially designed as ',...
            'a difficult upscaling test. Today, this highly heterogeneous ', ...
            'model is used a standard benchmark in a large number ',...
            'of publications. MRST has a dedicated module, "spe10", ', ...
            'that downloads and sets up the dataset'], ...
        'modelType', 'Two-phase fluid, Cartesian grid' ...
         );
     present = true;
end