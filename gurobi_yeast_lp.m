% Set path to Gurobi and load model
addpath(genpath('C:\gurobi910'))
load yeastGEM

% Set glucose and oxygen uptake rates
iGlu = find(strcmp(model.rxnNames,'D-glucose exchange'));
model.lb(iGlu) = -12.0;
iO2 = find(strcmp(model.rxnNames,'oxygen exchange'));
model.lb(iO2) = -8;

% Setup LP for Gurobi
LP.A          = sparse(model.S);
LP.obj        = model.c;
LP.rhs        = model.b;
LP.sense      = '='
LP.lb         = model.lb;
LP.ub         = model.ub;
LP.modelsense = 'max';

% Solve LP
result = gurobi(LP);

% Find growth and product secretion fluxes
iGr = find(strcmp(model.rxnNames,'growth'));
idxE = find(contains(model.rxnNames,'exchange'));
idxP = find(result.x(idxE) > 0);
idxFlux = [iGlu; iO2; iGr; idxE(idxP)];

% Create table to show LP results
table(model.rxnNames(idxFlux),result.x(idxFlux))