classdef OilWaterPolymerModel < TwoPhaseOilWaterModel
%
%
% SYNOPSIS:
%   model = OilWaterPolymerModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Two phase model with polymer. A description of the polymer model
% that is implemented here can be found in the directory ad-eor/docs .
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock structure
%   fluid    - Fluid structure
%   varargin - optional parameter
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: ThreePhaseBlackOilPolymerModel
%

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
    
    properties
        % Polymer present
        polymer
        
    end
    
    methods
        function model = OilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWPoly', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end
        
        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'Polymer', model.G.cells.num, 1);
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);
            
            if model.polymer
                c = model.getProp(state, 'polymer');
                c = min(c, model.fluid.cmax);
                state = model.setProp(state, 'polymer', max(c, 0) );
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'polymer');
                cmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cmax, c));
            end
        end

        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'polymer'}
                    index = 1;
                    fn = 'c';
                case {'polymermax'}
                    index = 1;
                    fn = 'cmax';
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                                    model, name);
            end
        end
        
        function scaling = getScalingFactorsCPR(model, problem, names)
            nNames = numel(names);

            scaling = cell(nNames, 1);
            handled = false(nNames, 1);

            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'polymer'
                        s = 0;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);

                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                % Get rest of scaling factors
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled));
                [scaling{~handled}] = other{:};
            end
        end

    end
end
