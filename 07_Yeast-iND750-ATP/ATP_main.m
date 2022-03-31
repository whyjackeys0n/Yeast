%% Load model
clear;
clc;
load iND750;
model = iND750;

%% find metabolite Sucrose
iScr = find(contains(model.rxnNames, 'Sucrose'));
model.rxnNames(iScr)

%% find reactions involved metabolite Sucrose
idx = find(model.S(930,:) ~= 0);
model.rxnNames(idx)

%% find metabolites involved transport reaction
idx1 = find(model.S(:,1174) ~= 0);
model.mets(idx1)
model.S(:,1174)
% sucr[e]+h[e] --> sucr[c]+h[c]
