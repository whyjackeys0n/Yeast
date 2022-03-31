function Gr = testoxygen(o2,strain)
% Set path to Gurobi and load model
addpath(genpath('C:\gurobi910'));

eval(['load ',strain]);
eval(['model = ',strain]);

% model.ub(1173) = 0;
% Set glucose and oxygen uptake rates
iGlu = find(strcmp(model.rxnNames,'D-Glucose exchange'));
model.lb(iGlu) = 0;
iFru = find(strcmp(model.rxnNames,'D-Fructose exchange'));
model.lb(iFru) = 0;
iScr = find(strcmp(model.rxnNames,'Sucrose exchange'));
model.lb(iScr) = -10;
iO2 = find(strcmp(model.rxnNames,'O2 exchange'));
model.lb(iO2) = -o2;
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
Gr = result.x(iGr);
end

