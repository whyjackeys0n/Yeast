%% Load model
clear;
clc;
load yeastGEM;

%% find metabolite Sucrose
iScr = find(contains(model.metNames, 'sucr'));
irxn = find(contains(model.rxnNames, 'sucr'));
model.metNames(iScr)

%% find reactions involved metabolite Sucrose
idx = find(model.S(1133,:) ~= 0);
model.rxnNames(idx)

%% find metabolites involved transport reaction
idx1 = find(model.S(:,783) ~= 0);
model.metNames(idx1)

% sucr[e]+h[e] --> sucr[c]+h[c]
