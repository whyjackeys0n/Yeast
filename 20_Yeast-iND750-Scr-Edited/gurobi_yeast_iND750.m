% Clear Workspace
clear;
clc;

% Set path to Gurobi and load model
addpath(genpath('C:\gurobi910'));
load iND750_Cheat;
model = iND750_Cheat;

% model.ub(1173) = 0;
% Set glucose and oxygen uptake rates
iGlu = find(strcmp(model.rxnNames,'D-Glucose exchange'));
model.lb(iGlu) = -10;
iFru = find(strcmp(model.rxnNames,'D-Fructose exchange'));
model.lb(iFru) = -10;
iScr = find(strcmp(model.rxnNames,'Sucrose exchange'));
model.lb(iScr) = 0;
iO2 = find(strcmp(model.rxnNames,'O2 exchange'));
model.lb(iO2) = -1;
iATP = find(strcmp(model.rxnNames,'ATP maintenance requirement'));
model.lb(iATP) = 0;

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
iGr = find(strcmp(model.rxnNames,'biomass SC4 bal'));
idxE = find(contains(model.rxnNames,'exchange'));
idxP = find(result.x(idxE) > 0);
idxFlux = [iGlu; iScr; iO2; iGr; idxE(idxP)];

% Create table to show LP results
table(model.rxnNames(idxFlux),result.x(idxFlux))

