function [problem, state] = equationsOilWaterDP(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for oil/water system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
G = model.G;
f = model.fluid;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0] = model.getProps(state0, 'pressure', 'water');

% Matrix properties
[pom,swm] = model.getProps(state, 'pom','swm');
[pom0,swm0] = model.getProps(state0, 'pom','swm');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, pom, swm, qWs, qOs, pBH] = initVariablesADI(p, sW, pom, swm, qWs, qOs, pBH);
    else
        zw = zeros(size(pBH));
        [p0, sW0, pom0, swm0, zw, zw, zw] = initVariablesADI(p0, sW0, pom0, swm0, zw, zw, zw); %#ok
        clear zw
    end
end
% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.
primaryVars = {'pressure', 'sW', 'pom', 'swm', 'qWs', 'qOs', 'bhp'};

%% Properties for Fracture
% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

%% Properties for Matrix
pvMultm = pvMult;
pvMultm0 = pvMult0;

% Using capillary pressure information
pcOWm = 0;
pcOWm0 = 0;
if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
    pcOWm  = model.fluid_matrix.pcOW(swm);
    pcOWm0  = model.fluid_matrix.pcOW(swm0);
end

pwm = pom - pcOWm;
pwm0 = pom0 - pcOWm0;

% SMALL TO DO HERE: WE USE THE REL PERMS OF THE MATRIX TO EVALUATE THE
% EFFECTIVE PERMEABILITY
som = 1-swm;
som0 = 1-swm0;
[krWm, krOm] = model.evaluteRelPerm({swm, som});

sigma = model.dp_info.shape_factor;
km = model.rock_matrix.perm(:,1);
muwm = model.fluid_matrix.muW(pwm);
muom = model.fluid_matrix.muO(pom);
bWm = f.bW(pwm);
bOm = f.bO(pom);

bWm0 = f.bW(pwm0);
bOm0 = f.bO(pom0);

%% Transfer
vb = model.G.cells.volumes;
Twm = model.dp_info.transfer_water(pW,pwm,sW,swm);
Tom = model.dp_info.transfer_oil(p,pom,sO,som);

%% Output Additional Info
if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, []);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
    state.Twm = double(Twm);
    state.Tom = double(Tom);
    state.pcm = double(pcOWm);
end

%% EQUATIONS ---------------------------------------------------------------
% Fracture
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.

bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Conservation of mass for water - fracture
water_fracture = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);
water_fracture = water_fracture + Twm;

% Conservation of mass for oil - fracture
oil_fracture = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);
oil_fracture = oil_fracture + Tom;

eqs{1} = water_fracture;
eqs{2} = oil_fracture;

% Add in any fluxes / source terms prescribed as boundary conditions.
eqs = addFluxesFromSourcesAndBC(model, eqs, ...
                                       {pW, p},...
                                       {rhoW,     rhoO},...
                                       {mobW,     mobO}, ...
                                       {bW, bO},  ...
                                       {sW, sO}, ...
                                       drivingForces);
                              
% Matrix 
% Conservation of mass for water - matrix
water_matrix = (s.pv/dt).*( pvMultm.*bWm.*swm - pvMultm0.*bWm0.*swm0 );
water_matrix = water_matrix - Twm;

% Conservation of mass for oil - fracture
oil_matrix = (s.pv/dt).*( pvMultm.*bOm.*som - pvMultm0.*bOm0.*som0 );
oil_matrix = oil_matrix - Tom;

eqs{3} = water_matrix;
eqs{4} = oil_matrix;

names = {'water_fracture', 'oil_fracture','water_matrix', 'oil_matrix'};
types = {'cell', 'cell', 'cell', 'cell'};                                   
                                   
%% Finally, add in and setup well equations
if ~isempty(W)
    wm = model.wellmodel;
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [model.fluid.rhoWS, model.fluid.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s = {sW(wc), sO(wc)};

        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration);
        % Store the well equations (relate well bottom hole pressures to
        % influx).
        eqs(5:6) = weqs;
        % Store the control equations (trivial equations ensuring that each
        % well will have values corresponding to the prescribed value)
        eqs{7} = ctrleqs;
        % Add source terms to the equations. Negative sign may be
        % surprising if one is used to source terms on the right hand side,
        % but this is the equations on residual form.
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};
%         
%         eqs{3}(wc) = eqs{3}(wc) - cqs{1};
%         eqs{4}(wc) = eqs{4}(wc) - cqs{2};
        
        names(5:7) = {'waterWells', 'oilWells', 'closureWells'};
        types(5:7) = {'perf', 'perf', 'well'};
    else
        [eqs(5:7), names(5:7), types(5:7)] = wm.createReverseModeWellEquations(model, state0.wellSol, p0);
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
