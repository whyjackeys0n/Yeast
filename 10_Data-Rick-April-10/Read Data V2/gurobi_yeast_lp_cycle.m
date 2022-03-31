clear
clc

% Set path to Gurobi and load model
load yeastGEM

Gr_exp = 0.34053;
Glu_start = 2;
Glu_end = 30;
for i = Glu_start:0.1:Glu_end
    % Set glucose and oxygen uptake rates
    iGlu = find(strcmp(model.rxnNames,'D-glucose exchange'));
    model.lb(iGlu) = -i;
    iO2 = find(strcmp(model.rxnNames,'oxygen exchange'));
    model.lb(iO2) = -8;

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
    idxFlux = [iGlu; iO2; iGr];
    
    % Find growth rate
    if result.x(iGr) > Gr_exp
        break;
    end
 
end
% Create table to show LP results
table(model.rxnNames(idxFlux),result.x(idxFlux))