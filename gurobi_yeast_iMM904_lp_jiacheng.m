% Clear Workspace
clear;
clc;

% Set path to Gurobi and load model
addpath(genpath('C:\gurobi910'));
load iMM904;
model = iMM904;

% model.ub(1173) = 0;
% Set glucose and oxygen uptake rates
iGlu = find(strcmp(model.rxnNames,'D-Glucose exchange'));
model.lb(iGlu) = 0;
iScr = find(strcmp(model.rxnNames,'Sucrose exchange'));
model.lb(iScr) = -10;
iO2 = find(strcmp(model.rxnNames,'O2 exchange'));
model.lb(iO2) = -1;
iATP = find(strcmp(model.rxnNames,'ATP maintenance requirement'));
model.lb(iATP) = 10;
model.ub(iATP) = 1000;
% iSuc = find(strcmp(model.rxnNames,'Sucrose exchange'));
% model.lb(iSuc) = -2;

% Setup LP for Gurobi
LP.A          = sparse(model.S);
LP.obj        = model.c;
LP.rhs        = model.b;
LP.sense      = '=';
LP.lb         = model.lb;
LP.ub         = model.ub;     
LP.modelsense = 'max';
 
% Solve LP
result = gurobi(LP);

% Find growth and product secretion fluxes
iGr = find(strcmp(model.rxnNames,'Biomass SC5 notrace'));
idxE = find(contains(model.rxnNames,'exchange'));
idxP = find(result.x(idxE) > 0);
idxFlux = [iScr; iO2; iGr; idxE(idxP)];

% Create table to show LP results
table(model.rxnNames(idxFlux),result.x(idxFlux))

% Find ethanol
% results = [];
% O2times = 5;
% iEthanol = find(strcmp(model.rxnNames,'Ethanol exchange'));
% for i = 1:O2times*(-model.lb(iGlu))
%     model.lb(iO2) = -i;
%     LP.lb = model.lb;
%     resulti = gurobi(LP);
%     resulti_Ethanol = resulti.x(iEthanol);
%     results = [results; resulti_Ethanol];
% end
% plot(1:O2times*(-model.lb(iGlu)), results);
