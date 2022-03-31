% Clear Workspace
clear;
clc;

% Set path to Gurobi and load model
addpath(genpath('C:\gurobi910'))
load yeastGEM

% Set glucose and oxygen uptake rates
iGlu = find(strcmp(model.rxnNames, 'D-glucose exchange'));
model.lb(iGlu) = -10;
iFru = find(strcmp(model.rxnNames, 'D-fructose exchange'));
model.lb(iFru) = -10;
iSuc = find(strcmp(model.rxnNames, 'sucrose exchange'));
model.lb(iSuc) = -10;
iO2 = find(strcmp(model.rxnNames, 'oxygen exchange'));
model.lb(iO2) = -1000;
iATP = find(strcmp(model.rxnNames, ...
    'non-growth associated maintenance reaction'));
model.lb(iATP) = 0.7;
model.ub(iATP) = 1000;

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
iGr = find(strcmp(model.rxnNames,'growth'));
idxE = find(contains(model.rxnNames,'exchange'));
idxP = find(result.x(idxE) > 0);
idxFlux = [iGlu; iFru; iSuc; iO2; iATP; iGr; idxE(idxP)];

% Create table to show LP results
table(model.rxnNames(idxFlux),result.x(idxFlux))

% Find ethanol
% results = [];
% O2times = 5;
% iEthanol = find(strcmp(model.rxnNames,'ethanol exchange'));
% for i = 1:O2times*(-model.lb(iGlu))
%     model.lb(iO2) = -i;
%     LP.lb = model.lb;
%     resulti = gurobi(LP);
%     resulti_Ethanol = resulti.x(iEthanol);
%     results = [results; resulti_Ethanol];
% end
% plot(1:O2times*(-model.lb(iGlu)), results);
