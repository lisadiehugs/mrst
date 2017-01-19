classdef TwoPhaseOilWaterDPModel < ThreePhaseBlackOilDPModel
% Two phase oil/water system without dissolution
properties

end

methods
    function model = TwoPhaseOilWaterDPModel(G, rock, fluid, rock_matrix, fluid_matrix, dp_info, varargin)
        model = model@ThreePhaseBlackOilDPModel(G, rock, fluid, rock_matrix, fluid_matrix, dp_info, varargin{:});

        % This is the model parameters for oil/water
        model.oil = true;
        model.gas = false;
        model.water = true;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        model.saturationVarNames = {'sw', 'so'};
        model.wellVarNames = {'qWs', 'qOs', 'bhp'};
        model.matrixVarNames = {'pom', 'swm'};
        
        model.extraStateOutput = 1;

        model = merge_options(model, varargin{:});
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsOilWaterDP(state0, state, model,...
                        dt, ...
                        drivingForces,...
                        varargin{:});

    end
end
end
