clear;
clc;
% Set path to Gurobi and load model
addpath(genpath('C:\gurobi910'))
load iND750
model = iND750;

% Set glucose and oxygen uptake rates
iGlu = find(strcmp(model.rxnNames,'D-Glucose exchange'));
model.lb(iGlu) = -19.4444*0.55;
iO2 = find(strcmp(model.rxnNames,'O2 exchange'));
model.lb(iO2) = -8*0.75;
iATP = find(strcmp(model.rxnNames,'ATP maintenance requirement'));
% model.lb(iATP) = 30;
% model.ub(iATP) = 1000;

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
idxFlux = [iGlu; iO2; iGr; idxE(idxP)];

% Create table to show LP results
table(model.rxnNames(idxFlux),result.x(idxFlux))