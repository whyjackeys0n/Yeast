clear;
clc;

load('iND750.mat');
model = iND750;

idx_scr_trans = find(contains(model.rxnNames,'sucrose transport'));
idx_scr_hydro = find(contains(model.rxnNames,'sucrose hydrolyzing'));

idx_ATP = find(strcmp(model.mets,'atp[c]'));
idx_ADP = find(strcmp(model.mets,'adp[c]'));
idx_Pi = find(strcmp(model.mets,'pi[c]'));
% idx_H2O = find(strcmp(model.mets,'h2o[c]'));
% idx_H = find(strcmp(model.mets,'h[c]'));

model.S([idx_ATP],idx_scr_trans) = -1;
model.S([idx_ADP,idx_Pi],idx_scr_trans) = 1;

idx_scr_hydro_h2oscr_e = find(model.S(:,idx_scr_hydro)<0);
idx_scr_hydro_fruglu_e = find(model.S(:,idx_scr_hydro)>0);
idx_scr_hydro_h2oscr_c = idx_scr_hydro_h2oscr_e - 1;
idx_scr_hydro_fruglu_c = idx_scr_hydro_fruglu_e - 1;
model.S(:,idx_scr_hydro) = 0;
model.S(idx_scr_hydro_h2oscr_c,idx_scr_hydro) = -1;
model.S(idx_scr_hydro_fruglu_c,idx_scr_hydro) = 1;

iND750_Private = model;
save iND750_Private iND750_Private

%% Cheat Model
model_cheat = iND750;
model_cheat.S(:,idx_scr_hydro) = 0;
model_cheat.S(:,idx_scr_trans) = 0;

iND750_Cheat = model_cheat;
save iND750_Cheat iND750_Cheat
